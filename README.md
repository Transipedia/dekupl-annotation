# DE-kupl Annotation

DE-kupl annotation is part of the DE-kupl package, and performs annotations of DE contigs identified by DE-kupl.

## Installation

Dependencies are cpan-minus (aka cpanm) and Dist::Zilla :

```
apt-get install cpanminus libdist-zilla-perl
```

The install rankvar with dzil and cpanm.

```
git clone https://gitlab.seq.one/workset/rankvar2.git && cd rankvar2
dzil install --install-command 'cpanm .'
```


## Dev environnement

For developpement, `git clone` this repository and enter the project folder.

Then, add the local dir to the PERL5LIB env var to use the modules locally.

```
export PERL5LIB=$PWD/lib:$PERL5LIB
```

You are ready to code!