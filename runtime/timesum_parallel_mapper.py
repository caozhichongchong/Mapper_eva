import os,glob
input_folder = '/media/caozhichongchong/baf7fe39-6567-421c-a79a-4c274b89f9f7/scripts/snp_curate/SNP_model_parallel_new/'
allinput_scripts = glob.glob('%s/*.out'%(input_folder))
allinput_scriptsmapper = glob.glob('%s/*mapper.sh.err'%(input_folder)) + glob.glob('%s/*mapper*mer.sh.err'%(input_folder))
allinput_scriptsmappermem = glob.glob('%s/*mappermem.sh.out'%(input_folder))

def get_thread(bashfile):
    thread_set = []
    tool = os.path.basename(bashfile).split('mapper')[0].split('.')[-1]
    for lines in open(bashfile):
        if '--num-threads' in lines:
            thread_set.append(tool + '\t' + lines.split('vcfer')[1].split(' ')[0] + '\t' + lines.split('--num-threads ')[1].split(' ')[0])
    return thread_set

def convert_to_seconds(time_str):
    # Split the string to separate minutes and seconds
    seconds = float(time_str.rstrip('s'))

    # Convert the time to seconds
    total_seconds = seconds
    return total_seconds

def convert_to_seconds_usertime(time_str):
    # Split the string to separate minutes and seconds
    parts = time_str.split(':')
    minutes = int(parts[-2])
    seconds = float(parts[-1])

    # Convert the time to seconds
    total_seconds = minutes * 60 + seconds
    return total_seconds

alltimeindex = ['genome\ttool\tthread\ttime\n']
alltimemapping = ['genome\ttool\tthread\ttime\n']
allmem = ['genome\ttool\tthread\tmemory\n']

for input_script in allinput_scripts:
        scriptname = os.path.basename(input_script)
        genome = scriptname.split('.')[0]
        tool = scriptname.split('.')[1]
        thread_set = [1,3,5,7,9,11,15,20,25,30]
        time_set_index = []
        time_set_mapping = []
        mem_set_index = []
        mem_set_mapping = []
        i = 0
        for lines in open(input_script, 'r'):
            if 'mapper' not in input_script:
                if ('Elapsed (wall clock) time (h:mm:ss or m:ss): ') in lines:
                    realtime = (lines.split('Elapsed (wall clock) time (h:mm:ss or m:ss): ')[1].split('\n')[0].replace(' ',''))
                    if i%2 == 0:
                        time_set_index.append(convert_to_seconds_usertime(realtime))
                    else:
                        if 'strobealign' in input_script:
                            # index cannot be used
                            time_set_mapping.append(convert_to_seconds_usertime(realtime) - time_set_index[-1])
                        else:
                            time_set_mapping.append(convert_to_seconds_usertime(realtime))
            if ('Maximum resident set size (kbytes): ') in lines:
                realmem = (lines.split('Maximum resident set size (kbytes): ')[1].split('\n')[0].replace(' ',''))
                if 'mapper' not in input_script:
                    if i%2 == 0:
                        mem_set_index.append(realmem)
                    else:
                        mem_set_mapping.append(realmem)
                elif 'mappermem' not in input_script:
                    print(realmem)
                    mem_set_index.append(realmem)
                    mem_set_mapping.append(realmem)
                i += 1
        if 'mapper' not in input_script:
            if 'last' in input_script:
                print(genome, tool, 'last convert maf to sam takes',time_set_index[-1], 'seconds')
                time_set_index = time_set_index[:-1]
            alltimeindex += ['%s\t%s\t'%(genome,tool) +'%s\t%s\n'%(x,y) for x,y in zip(thread_set,time_set_index)]
            alltimemapping += ['%s\t%s\t' % (genome, tool) + '%s\t%s\n' % (x, y) for x, y in zip(thread_set, time_set_mapping)]
        mem_set_all = [max(int(x),int(y))/1000000 for x,y in zip(mem_set_index, mem_set_mapping)]
        allmem += ['%s\t%s\t' % (genome, tool) + '%s\t%s\n' % (x, y) for x, y in zip(thread_set, mem_set_all)]
for input_script in allinput_scriptsmapper:
    scriptname = os.path.basename(input_script)
    genome = scriptname.split('.')[0]
    tool = scriptname.split('.')[1]
    thread_set = [1, 3, 5, 7, 9, 11, 15, 20, 25, 30]
    time_set_index = []
    time_set_mapping = []
    i = 0
    for lines in open(input_script, 'r'):
        if lines.startswith('Hashed reference '):
            if ' 54 ' in lines or ' 58 ' in lines:
                realtime = (lines.split('at ')[1].split('\n')[0].replace(' ', ''))
                time_set_index.append(convert_to_seconds(realtime))
        if lines.startswith('Done in '):
            realtime = (lines.split('Done in ')[1].split('.\n')[0].replace(' ', ''))
            print('mapper',convert_to_seconds(realtime),time_set_index[-1])
            time_set_mapping.append(convert_to_seconds(realtime)-time_set_index[-1])
            i += 1
    print('mapper',time_set_index,time_set_mapping)
    alltimeindex += ['%s\t%s\t' % (genome, tool) + '%s\t%s\n' % (x, y) for x, y in zip(thread_set, time_set_index)]
    alltimemapping += ['%s\t%s\t' % (genome, tool) + '%s\t%s\n' % (x, y) for x, y in zip(thread_set, time_set_mapping)]

f1 = open('%s/alltimeindex.txt'%(input_folder),'w')
f1.write(''.join(alltimeindex))
f1.close()
f1 = open('%s/alltimemapping.txt'%(input_folder),'w')
f1.write(''.join(alltimemapping))
f1.close()
f1 = open('%s/allmem.txt'%(input_folder),'w')
f1.write(''.join(allmem))
f1.close()

alltimemem = ['genome\ttool\tthread\tmemory_limit\ttime\tmemory\n']
for input_script in allinput_scriptsmappermem:
        scriptname = os.path.basename(input_script)
        genome = scriptname.split('.')[0]
        tool = scriptname.split('.')[1]
        thread_set = [5, 30]
        mem_set = [5,10,15,20,30]
        thread_mem_set = []
        for i in thread_set:
            for mem in mem_set:
                thread_mem_set.append('%s\t%s'%(i,mem))
        time_set = []
        mem_set = []
        for lines in open(input_script, 'r'):
            if ('Elapsed (wall clock) time (h:mm:ss or m:ss): ') in lines:
                realtime = (lines.split('Elapsed (wall clock) time (h:mm:ss or m:ss): ')[1].split('\n')[0].replace(' ',''))
                time_set.append(convert_to_seconds_usertime(realtime))
            if ('Maximum resident set size (kbytes): ') in lines:
                realmem = (lines.split('Maximum resident set size (kbytes): ')[1].split('\n')[0].replace(' ',''))
                mem_set.append(realmem)
        if 'last' in input_script:
            print(genome, tool, 'last convert maf to sam takes', mem_set[:-1], 'kbytes')
            mem_set = mem_set[:-1]
        mem_set_all = ['%s\t%s'%(x,y) for x,y in zip(time_set, mem_set)]
        alltimemem += ['%s\t%s\t' % (genome, tool) + '%s\t%s\n' % (x, y) for x, y in zip(thread_mem_set, mem_set_all)]
f1 = open('%s/allmappermemtime.txt'%(input_folder),'w')
f1.write(''.join(alltimemem))
f1.close()