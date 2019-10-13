import datetime
import paramiko
import sys
import threading
import time
from scp import SCPClient

############### configurations starts ###################
# instance config
pem_path = '/home/bn7/Documents/pem_files/WGCNA_Instance_KP.pem'
instance_dns = ['ec2-18-222-233-106.us-east-2.compute.amazonaws.com', 'ec2-3-16-22-102.us-east-2.compute.amazonaws.com', 'ec2-18-218-1-94.us-east-2.compute.amazonaws.com', 'ec2-3-16-43-173.us-east-2.compute.amazonaws.com', 'ec2-18-220-141-109.us-east-2.compute.amazonaws.com']
username = 'ubuntu'

# credentials config
# bitbucket
bb_uname = 'kcbg475'
bb_pwd = ''

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
def run_job(instance_dns, username, pem_path, id_path):
    fileprint("running job "+id_path)
    fileprint("instance: "+instance_dns)
    fileprint(id_path)
    client = paramiko.SSHClient()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    client.connect(
                    instance_dns,
                    username=username,
                    key_filename=pem_path)
    # send file with ids to the instance
    with SCPClient(client.get_transport()) as scp:
        scp.put(id_path, '~/'+id_path)
    # clone the repo in the instance using git
    stdin, stdout, stderr = client.exec_command('git clone https://'+bb_uname+':'+bb_pwd+'@bitbucket.org/elucidatainc/masswgcna.git')
    exit_status = stdout.channel.recv_exit_status()          # Blocking call
    # upgrade apt and install pip & aws 
    # install r and it's packages too
    stdin, stdout, stderr = client.exec_command('cd masswgcna; ./configure.sh > output-configure.out')
    exit_status = stdout.channel.recv_exit_status()          # Blocking call
    # configure aws key
    stdin, stdout, stderr = client.exec_command('aws configure set aws_access_key_id '+ my_access_key)
    stdin, stdout, stderr = client.exec_command('aws configure set aws_secret_access_key '+ my_secret_key)
    # run r script here
    stdin, stdout, stderr = client.exec_command('cd masswgcna; ./masswgcna.sh -t '+ dataset_type +' -g ~/'+id_path+' -f '+output_folder+' > output-masswgcna.out')
    exit_status = stdout.channel.recv_exit_status()          # Blocking call
    # close connection
    client.close()

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
