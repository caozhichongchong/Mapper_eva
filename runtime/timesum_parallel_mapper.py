import os,glob
input_folder = '/media/caozhichongchong/baf7fe39-6567-421c-a79a-4c274b89f9f7/scripts/snp_curate/SNP_model_parallel_new/'
allinput_scripts = glob.glob('%s/*.out'%(input_folder))
def get_thread(bashfile):
    thread_set = []
    tool = os.path.basename(bashfile).split('mapper')[0].split('.')[-1]
    for lines in open(bashfile):
        if '--num-threads' in lines:
            thread_set.append(tool + '\t' + lines.split('vcfer')[1].split(' ')[0] + '\t' + lines.split('--num-threads ')[1].split(' ')[0])
    return thread_set

def convert_to_seconds(time_str):
    # Split the string to separate minutes and seconds
    parts = time_str.split('m')
    minutes = int(parts[0])
    seconds = float(parts[1].rstrip('s'))

    # Convert the time to seconds
    total_seconds = minutes * 60 + seconds
    return '%s\t%s'%(time_str,total_seconds)

alltime = ['genome\ttool\tthread\treal time\ttime\n']
for input_script in allinput_scripts:
    scriptname = os.path.basename(input_script)
    genome = scriptname.split('.')[0]
    tool = scriptname.split('.')[1]
    thread_set = range(1, 30, 2)
    time_set = []
    for lines in open(input_script, 'r'):
        if lines.startswith('real'):
            realtime = (lines.split('real')[1].split('\n')[0].replace(' ',''))
            time_set.append(convert_to_seconds(realtime))
    alltime += ['%s\t%s\t'%(genome,tool) +'%s%s\n'%(x,y) for x,y in zip(thread_set,time_set)]
f1 = open('%s/alltime.txt'%(input_folder),'w')
f1.write(''.join(alltime))
f1.close()