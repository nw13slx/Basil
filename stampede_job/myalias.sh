#my alias setting on Stampede, it is sourced when logging in
alias hw='echo hello word'

#job related
wmj() { 
    showq |grep "Running"|grep $MYNAME
    showq -l|grep "normal"|grep "Waiting"|cat -b |grep $MYNAME
    showq  -l|grep "normal"|grep "Waiting"|tail -n 10
}
alias sdev='showq -l|grep development'
alias sser='showq -l|grep serial'
alias curious='showq -l |grep "Running"|sort -k 5 -n'
alias mj='showq -l -u'


#rsync is not suggested for ranch.
#alias backup_rsync='rsync -avzh $SCRATCH $MYNAME@$ARCHIVER:$ARCHIVE'

alias fractal='ssh -X $MYNAME@fractal.mit.edu'
alias ranch='ssh $MYNAME@ranch.tacc.utexas.edu'
alias sftphopper='sftp $MYNAME@hopper.nersc.gov'
alias sftpedison='sftp $MYNAME@edison.nersc.gov'
alias sftpcori='sftp $MYNAME@cori.nersc.gov'
alias sftpfractal='sftp $MYNAME@fractal.mit.edu'
alias sftpranch='sftp $MYNAME@ranch.tacc.utexas.edu'

alias sls='ls -t -n|sort -k 5 -n'
alias lst='ls -Srtnlh'

#VASP Related
alias cleanvasp='rm LOCPOT CHG* EIGENVAL OSZICAR WAVECAR DOSCAR IBZKPT PCDAT REPORT vasprun.xml XDATCAR PROCAR'
alias dos='~/bin/split_dos.py $(awk "NR==7" CONTCAR |wc -w) $(awk "NR==7" CONTCAR)'
alias split='~/bin/spin_density CHGCAR'
baderana () {
    ~/bin/spin_density CHGCAR
    mkdir elec mag
    cd elec
    ~/bin/bader ../CHGCAR
    cd ../mag
    ~/bin/bader ../CHGCAR-u-d -ref ../CHGCAR
    cd ../
    rm CHGCAR-u-d
}
baderanavac () {
    ~/bin/spin_density CHGCAR
    mkdir elec mag
    cd elec
    ~/bin/bader -vac auto ../CHGCAR
    cd ../mag
    ~/bin/bader -vac auto ../CHGCAR-u-d -ref ../CHGCAR
    cd ../
    rm CHGCAR-u-d
}
