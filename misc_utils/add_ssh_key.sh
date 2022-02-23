ssh-keygen -t rsa
ssh-copy-id -o ProxyCommand="ssh -W %h:%p blasseln@cermics.enpc.fr" blasseln@clustern14
exec ssh-agent sh
