FROM registry.gitlab.com/transipedia/dekupl-annotation:base

WORKDIR /dekupl
COPY . .
RUN dzil install --install-command 'cpanm .'

ENTRYPOINT [ "dkpl" ]
