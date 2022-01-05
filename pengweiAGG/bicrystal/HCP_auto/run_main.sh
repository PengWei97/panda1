source /home/opt/anaconda3/bin/activate
# conda activate moose
nohup python -u main.py > output.out 2>&1 & echo $! > run.pid