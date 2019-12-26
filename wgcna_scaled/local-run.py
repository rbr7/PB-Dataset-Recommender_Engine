import datetime
import paramiko
import sys
import threading
import time
from scp import SCPClient

############### configurations starts ###################
# instance config
pem_path = ##
instance_dns = ##
username = ##

# credentials config
# ##
bb_uname = ##
bb_pwd = ##

# aws
my_access_key = ''
my_secret_key = ''
############### configurations ends ###################

# folder/file/run specifications
master_id_file = sys.argv[1]
output_folder = 'master_run_'+str(datetime.datetime.now())
dataset_type = sys.argv[2]

def fileprint(message):
    with open("output.txt", 'a') as fp:
        fp.write(message+'\n')
        fp.close()

# this is the job to be run on each instance
###################################

# function to divide major file into x smaller files where x is the number of instances
def divide_id_file(id_file, div_len=len(instance_dns)):
    id_paths = []
    with open(id_file, 'r') as fp:
        lines = fp.readlines()
        each_len =int( len(lines)/div_len)
        for i in range(0,div_len):
            each = lines[i*each_len: (i+1)*each_len]
            with open('each_'+str(i), 'w') as wfp:
                wfp.writelines(each)
                wfp.close()
                id_paths.append('each_'+str(i))
        fp.close()
    return id_paths


# main function
id_paths = divide_id_file(master_id_file)
print("sending jobs to instances")
jobs = []
threads = []
for instance, id_path in zip(instance_dns, id_paths):
    try:
        print(id_path)
        thread = threading.Thread(target=run_job, args=(instance, username, pem_path, id_path), 
        daemon=False)
        threads.append(thread)
        thread.start()
    except Exception as e:
        print(e)

for thread in threads:
    thread.join()
