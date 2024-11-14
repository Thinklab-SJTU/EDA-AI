user=$(whoami)
ps -o pid,etime,cmd -u $user | grep main.py | grep python | grep circuit > ps.log

# print each line in ps.log, after output each line, add an \n
num_process=0
while read line
do
    echo $line
    echo
    num_process=$(($num_process+1))
done < ps.log

echo "Number of processes: $num_process"

if [ -f ps.log ]
then
    rm ps.log
fi
