package DEkupl::Analyzer::Switches;
# ABSTRACT: Compute Differential usage between contig and gene

use Moose;

use DEkupl::Utils;
use File::Temp;

with 'DEkupl::Analyzer';

has 'sample_names' => (
  is => 'ro',
  isa => 'ArrayRef[Str]',
  required => 1,
  traits => ['Array'],
  handles => {
    all_samples => 'elements',
  },
);

has 'normalized_gene_counts_file' => (
  is => 'ro',
  isa => 'Str',
  required => 1,
);

has 'sample_conditions_file' => (
  is => 'ro',
  isa => 'Str',
  required => 1,
);

# TODO Verify that input files are valid!

my @columns = ('du_pvalue', 'du_stat');

sub BUILD {
  my $self = shift;

  # TODO Verify that R is available as well as the packages we need.
  # TODO Split this into functions and test them!
  my $check_log = new File::Temp(SUFFIX => '.log', UNLINK => 0);
  if(system("R --slave -e 'library(DESeq2);' 2> $check_log") != 0) {
    die "R and/or DESeq2 package are not available (see logs in $check_log)";
  } else {
    unlink $check_log;
  }

  # Verify that the sample condition file is correct
  verifySampleConditionsFile($self->sample_conditions_file, $self->sample_names);

  # Create a temp file with the R script
  my $rscript = new File::Temp( SUFFIX => '.R', UNLINK => 0);

  print $rscript getRscript();
  close $rscript;

  $self->verboseLog("Temporary Rscript file created at $rscript");

  # Change permission to execute the script
  #chmod 0755, $rscript;

  # Now we create a temp file with samples counts per contigs and gene information
  my $contigs_file = new File::Temp( SUFFIX => '.tsv', UNLINK => 0);

  # Print headers
  print $contigs_file join("\t",
        "contig_id",
        "gene_id",
        @{$self->sample_names},
  ), "\n";

  # Store gene_ids into memory
  my %gene_ids;

  my $contigs_it = $self->contigs_db->contigsIterator();
  my $nb_contigs_with_genes = 0;
  while(my $contig = $contigs_it->()) {
    my $gene_id = $contig->{gene_id};
    if(defined $gene_id) {
      $nb_contigs_with_genes++;
      $gene_ids{$gene_id} = 1;
      print $contigs_file join("\t",
        $contig->{tag},
        $contig->{gene_id},
        @{$contig->{counts}},
      ), "\n";
    }
  }
  $self->verboseLog("$nb_contigs_with_genes contigs with genes will be scanned for switches");

  $self->verboseLog("Temporary contigs file created at $contigs_file");

  close $contigs_file;

  # Now we create a temp file with samples counts per contigs and gene information
  my $genes_file = new File::Temp( SUFFIX => '.tsv', UNLINK => 0);

  {
    # Print headers
    print $genes_file join("\t",
        "gene_id",
        $self->all_samples
    ), "\n";

    my $gene_counts_fh = DEkupl::Utils::getReadingFileHandle($self->normalized_gene_counts_file);
    my $header_line = <$gene_counts_fh>;
    chomp $header_line;
    my @headers = split "\t", $header_line;
    $headers[0] = "gene_id"; # Force the name of the feature column

    my $nb_genes_with_contigs = 0;

    while(<$gene_counts_fh>) {
      chomp;
      my @values = split "\t", $_;
      my %fields = map { $headers[$_] => $values[$_] } 0..$#headers;

      my $gene_id = DEkupl::Utils::getAtomicGeneID($fields{gene_id});

      #print STDERR "GENE_ID $gene_id\n"; sleep 1;

      # Only print genes that have DE contigs
      next if !defined $gene_ids{$gene_id};

      $nb_genes_with_contigs++;

      print $genes_file join("\t",
        $gene_id,
        map { $fields{$_} } $self->all_samples
      ), "\n";
    }

    $self->verboseLog("$nb_genes_with_contigs genes with associated contigs found in DEGs");
  }

  $self->verboseLog("Temporary genes file created at $genes_file");

  close $genes_file;

  # We execute the getSwitches script
  my $switches_file = new File::Temp( SUFFIX => '.tsv', UNLINK => 0);
  close($switches_file);

  # Temp file to place logs from R execution
  my $rscript_logs = new File::Temp(SUFFIX => '.log', UNLINK => 0);

  my $command = join(" ",
    "Rscript",
    $rscript,
    $switches_file,
    $contigs_file,
    $self->sample_conditions_file,
    $genes_file,
    "2>",
    $rscript_logs
  );

  $self->verboseLog("Command: $command");

  system($command) == 0 or die ("Error(s) in Switch computing using R and DESeq2 (see logs in $rscript_logs)");

  $self->verboseLog("Temporary switches file created at $switches_file");

  # Everything went well, we remove the logs and input files
  unlink $rscript_logs;
  unlink $rscript;
  unlink $contigs_file;
  unlink $genes_file;

  {
    my $switches_fh = DEkupl::Utils::getReadingFileHandle($switches_file);
    my $header_line = <$switches_fh>;
    while(<$switches_fh>){
      chomp;
      my ($contig_id, $du_pvalue, $du_stat) = split "\t", $_;
      my $contig = $self->contigs_db->loadContig($contig_id);
      $contig->{du_pvalue} = $du_pvalue;
      $contig->{du_stat} = $du_stat;
      # Save contig
      $self->contigs_db->saveContig($contig);
    }
  }

}

sub verifySampleConditionsFile {
  my $file          = shift;
  my $sample_names  = shift;
  my $fh = DEkupl::Utils::getReadingFileHandle($file);
  my $header_line = <$fh>;
  my %samples_check = map { $_ => 0 } @{$sample_names};
  my $i = 1;
  while(<$fh>) {
    chomp;
    my ($sample, $condition, $normalization_factor) = split "\t", $_;
    die("Line $i did not had three columns in file $file (file must be tabulated") unless defined $sample && defined $condition && defined $normalization_factor;
    die("Sample name $sample does not match contig file at line $i in file $file") unless defined $samples_check{$sample};
    die("Normalization factor ($normalization_factor) must have a positive value at line $i in file $file") unless $normalization_factor > 0;
    $samples_check{$sample}++;
    $i++;
  }
  # Verify is all samples in the contigs file have been found and that there is no duplicated entries
  foreach my $k (keys %samples_check) {
    if($samples_check{$k} == 0) {
      die("Sample $k not found in file $file");
    } elsif($samples_check{$k} > 1) {
      die("Duplicated entry for sample $k in file $file");
    }
  }
  return 0;
}

sub getHeaders {
  my $self = shift;
  return @columns;
}

sub getValues {
  my $self = shift;
  my $contig = shift;
  my @values = map { defined $contig->{$_}? $contig->{$_} : 'NA' } @columns;
  return @values;
}

sub getRscript {
  my $rscript_content =<<END;
#!/usr/bin/env Rscript

args <- commandArgs(TRUE)

if (length(args)==0){
  stop("missing arguments !\nUsage : ./getSwitches.R <output_file.tsv> <contigs.tsv> <sample_conditions.tsv> <normalized_counts.tsv>")
}

library(DESeq2)

# Get files for args

output_file                 <- args[1] # Output file
all_contigs_file            <- args[2] # DiffContigsInfos from dekupl
sample_conditions_file      <- args[3] # design file from dekupl
normalized_gene_counts_file <- args[4] # normalized counts (gene expression) from dekupl (Kallisto)

# Open files
tab_counts_DEkupl     <- read.delim(all_contigs_file,check.names=F)
sample_conditions     <- read.delim(sample_conditions_file,check.names=F)
normalizedGeneCounts  <- read.delim(normalized_gene_counts_file,check.names=F)

#### process data ####

#intersect KALLISTO gene IDs & DEKUPL gene IDs
tab_counts_DEkupl     <- tab_counts_DEkupl[which(tab_counts_DEkupl\$gene_id \%in% normalizedGeneCounts\$gene_id),]
normalizedGeneCounts  <- normalizedGeneCounts[which(normalizedGeneCounts\$gene_id \%in% tab_counts_DEkupl\$gene_id),]
tab_counts_Kallisto   <- merge(tab_counts_DEkupl,normalizedGeneCounts, by.x="gene_id", by.y="gene_id", all.x=T, all.y=F)

#reorganize columns in order to have contig ID, gene ID & KALLISTO counts (same order as tab_counts_DEkupl)
tab_counts_Kallisto   <- tab_counts_Kallisto[,c(2,1,(ncol(tab_counts_DEkupl)+1):(ncol(tab_counts_Kallisto)))]

#order both tables following the contig ID
tab_counts_Kallisto   <- tab_counts_Kallisto[order(tab_counts_Kallisto\$contig_id),]

tab_counts_DEkupl     <- tab_counts_DEkupl[order(tab_counts_DEkupl\$contig_id),]

#keep the same header for both tables
names(tab_counts_Kallisto)[3:ncol(tab_counts_Kallisto)] <- names(tab_counts_DEkupl)[3:ncol(tab_counts_DEkupl)]

#prepare contigs with their counts for DESeq2 (row names = contig ID, and we keep only counts without any other columns)
rownames(tab_counts_DEkupl) <- tab_counts_DEkupl\$contig_id

tab_counts_DEkupl[,c(1,2)] <- NULL

#get conditions name
cond1 <- as.character(sample_conditions[1,2])
cond2 <- unique(as.character(sample_conditions[,2][sample_conditions[,2]!=cond1]))

#get number of samples for each condition
rep_cond1 <- nrow(sample_conditions[which(sample_conditions[,2]==cond1),])

rep_cond2 <- nrow(sample_conditions[which(sample_conditions[,2]==cond2),])

#set design
samples <-  data.frame(row.names=names(tab_counts_DEkupl),condition=c(rep(cond1,rep_cond1),rep(cond2,rep_cond2)))

#create DESeqDataSet object from design & contigs DE
dds <- DESeqDataSetFromMatrix(countData=as.matrix(round(tab_counts_DEkupl)),colData=samples,design=~condition)

#compute normalization factor for each contig at each sample, thanks to their normalized gene counts from Kallisto
normFactors <- as.matrix((tab_counts_Kallisto[,3:ncol(tab_counts_Kallisto)]+1)/exp(rowMeans(log(tab_counts_Kallisto[,3:ncol(tab_counts_Kallisto)]+1))))

#allocation of normalization factors
normalizationFactors(dds) <- normFactors

#estimate overdispersion parameters
#it's possible to have issues with estimateDispersions() if you have a low number of contigs ("dispersion trend not well captured")
#so, we use fitType="mean" instead of the default "parametric"
getDispersions <- function(my_object="") {

  dds<-try(estimateDispersions(my_object))

  if (class(dds)=="try-error"){
    cat("with fitType='parametric', the dispersion trend was not well captured by the function, we will use instead fitType='mean'")
    dds<-estimateDispersions(my_object,fitType="mean")
  }

  return(dds)
}

dds <- getDispersions(dds)
dds <- nbinomWaldTest(dds) #binomiale negative test

#results
#we turn off all kind of filters to avoid "NA" values for outliers
DESeq2Result <- results(dds,independentFiltering=F,cooksCutoff=F)
DESeq2Result <- DESeq2Result[c("padj","stat")] # extract padj

#make a custom table with contig ID,mean cond1, mean cond2, log2FC, padj, normalized counts for all libraries
new_result <- data.frame(
  contig_id = row.names(DESeq2Result),
  du_pvalue = DESeq2Result\$padj,
  du_stat   = DESeq2Result\$stat,row.names=NULL
)

# Write the table to the output file
write.table(new_result, file=output_file, sep="\t", row.names=F, col.names=T, quote=F)
END
  return $rscript_content;
}

no Moose;
__PACKAGE__->meta->make_immutable;
