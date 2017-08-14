#my alias setting on Stampede, it is sourced when logging in
alias hw='echo hello word'

#job related
wmj() { 
    showq |grep "Running"|grep $USER
    showq -l|grep "normal"|grep "Waiting"|cat -b |grep $USER
    showq  -l|grep "normal"|grep "Waiting"|tail -n 10
}
alias sdev='showq -l|grep development'
alias sser='showq -l|grep serial'
alias curious='showq -l |grep "Running"|sort -k 5 -n'
alias mj='showq -l -u'
alias touchall="find . -exec touch {} \;"



#rsync is not suggested for ranch.
#alias backup_rsync='rsync -avzh $SCRATCH $USER@$ARCHIVER:$ARCHIVE'

alias fractal='ssh -X $USER@fractal.mit.edu'
alias ranch='ssh $USER@ranch.tacc.utexas.edu'
alias sftphopper='sftp $USER@hopper.nersc.gov'
alias sftpedison='sftp $USER@edison.nersc.gov'
alias sftpcori='sftp $USER@cori.nersc.gov'
alias sftpfractal='sftp $USER@fractal.mit.edu'
alias sftpranch='sftp $USER@ranch.tacc.utexas.edu'
alias sftpmghpcc='sftp $USER@eofe7.mit.edu'

alias sls='ls -t -n|sort -k 5 -n'
alias lst='ls -Srtnlh'

#VASP Related
alias dos='~/bin/split_dos.py $(awk "NR==7" CONTCAR |wc -w) $(awk "NR==7" CONTCAR)'
alias split='~/bin/spin_density CHGCAR'
