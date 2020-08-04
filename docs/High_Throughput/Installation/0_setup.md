# Getting started with the high-throughput version of FHI-vibes
FireWorks is a high-throughput framework that can greatly enhance the efficiency of your work, but does require initial setup to get working.
Because of this, it is important to consider if running calculations in high-throughput is necessary for your work, or if using the command line interface tools vibes provides is sufficient.
If you only need to do an in-depth study of a few materials then high-throughput aspects of FHI-vibes is probably not necessary, but if you need to do a systematic study over many materials this could be a useful tool for your research.

This setup guide will only focus on the basic aspects that you need do in order to run our FireWorks workflows, for more advanced usage please look at the [FireWorks online documentation](https://materialsproject.github.io/fireworks/).


## Setting up MongoDB
Before running or installing FireWorks you first need to have access to a MongoDB instance to act as the job management database (LaunchPad in FireWorks terminology).
To do this there are two options: set up your own database by installing MongoDB locally or on an accessible server, or by using a cloud provider.
If you want to install MongoDB yourself please follow [their instructions](https://www.mongodb.com/), for a cloud provider you can try [mlab](http://mlab.com/), [MongoDB Atlas](https://www.mongodb.com/cloud/atlas), or another option.

Once you have created a database you'll need to set up a user to access it.
Furthermore it is recommended that you also set up secure passwords for all your databases for data security reasons.
To set this up we suggest you follow this guide, but other methods are possible.
We have a separate admin/user structure in order to ensure that admin access to the database can't be obtained through a FireWorks configuration file.

### How to create a new MongoDB with authorization
*Note: This is for a locally managed database, if you are using a cloud service security procedures maybe different please consult their instructions*

First create a database with the correct base directory (\$BASE_DIR), port (\$PORT), and binding IP addresses (\$BIND_IP) with the following command
```
mongod --logpath $BASE_DIR/logs --dbpath $BASE_DIR/db --port $PORT --bind_ip $BIND_IP --fork
```
It is good practice to include 127.0.0.1 (localhost) as one of the binding IP's along with the one you will use to access the database from the outside.
The list of IP addresses can be stored in $BIND_IP with the following command
```
export BIND_IP=127.0.0.1,IP_1,IP_2,...
```
note it is a comma separated list.

From here you can access the database with
```
mongo --port $PORT --host $HOST_NAME
```
\$HOST_NAME should be one of the IP addresses in \$BIND_IP, if localhost is on the list and you are on that server, the host keyword is not necessary.
You should now be in a MongoDB terminal, from here type in the following commands
```
use admin
db.createUser(
  {
    user: "admin_user",
    pwd: "ADMIN_PASSWORD",
    roles: [ { role: "userAdminAnyDatabase", db: "admin" }, "readWriteAnyDatabase" ]
  }
)
exit
```
This will create an admin user that has global access to all databases on the server. Most importantly now you can restart your database requiring authorization to access it.
```
mongod --dbpath $BASE_DIR/db --shutdown
mongod --logpath $BASE_DIR/logs --dbpath $BASE_DIR/db --port $PORT --bind_ip $BIND_IP --fork --auth
```
Now that the database is running you can create additional users/databases inside the main one with the following commands.
To start launch the MongoDB terminal
```
mongo --port $PORT --host $HOST_NAME -u "admin_user" -p "ADMIN_PASSORD" --authenticationDatabase "admin"
```
then inside that terminal type the following commands
```
use FIREWORKS_DB_NAME
db.createUser(
  {
    user: "USER",
    pwd: "PASSWORD",
    roles: [ { role: "dbOwner", db: "FIREWORKS_DB_NAME" } ]
  }
)
exit
```
Please do not use the same username/password for the admin user and individual users.
Now that you have a database running you can store this information in a yaml file called my_launchpad.yaml with the following contents
```
host: HOST
port: PORT
name: FIREWORKS_DB_NAME
username: USER
password: PASSWORD
```
The host this time should be the IP address in $BIND_IP that you intend to use for external connections. If this is on the same machine you will be running the workflows you maybe able to use localhost, but check with your sysadmin on how the nodes are set up.
Now that you have a database set up now it is time to install FireWorks

## Installing FireWorks
The first step in using FireWorks is installing the python library.
By installing vibes with the FireWorks dependency by typing `poetry install -E fireworks` or `pip install vibes[fireworks]` this is already included, but you may want to install your own version of FireWorks locally.
You can do this either by cloning the [git repository](https://github.com/materialsproject/fireworks) and using the python setup tools or via pip/conda.
In addition to FireWorks, if you want to use the remote clients/database you'll need to install paramiko and fabric or if you are using the NEWT queuing system you'll have to install requests.
To do all of this with pip simply type in
```
pip install FireWorks
pip install paramiko  # (only needed if using built-in remote file transfer!)
pip install fabric  # (only needed if using daemon mode of qlaunch!)
pip install requests  # (only needed if you want to use the NEWT queue adapter!)
```
The fireworks extension of vibes already includes paramiko and fabric, but requests would have to be installed via pip.

Once FireWorks is installed you can test your installation by attempting to connect to FireWorks read-only test database, by creating `my_launchpad_testing.yaml` with the following contents:
```
host: ds049170.mongolab.com
port: 49170
name: fireworks
username: test_user
password: testing123
```
Then you can run the following command

```
lpad -l my_launchpad_testing.yaml get_wflows
```
And you should get the following output
```
[
    {
        "name": "Tracker FW--1",
        "state": "READY",
        "states_list": "REA",
        "created_on": "2014-10-27T15:00:25.408000"
    },
    {
        "name": "Tracker FW--2",
        "state": "READY",
        "states_list": "REA",
        "created_on": "2014-10-27T15:00:25.775000"
    }
]
```

Now that FireWorks is installed properly it is time to set up your configuration.
Create a fireworks configuration directory (\$FW_CONFIG) to store all the configuration files.
We recommend you use `.fireworks/` for \$FW_CONFIG.
Move the my_launchpad.yaml file to \$FW_CONFIG.
If you are planning to use FireWorks with a queuing system also create a my_qadapter.yaml file in \$FW_CONFIG. To get a sample of what to do for your queue system go to the [FireWork's git hub](https://github.com/materialsproject/fireworks/tree/master/fw_tutorials/queue) page and download the correct one.
For example here is one for a SLURM System
```
_fw_name: CommonAdapter
_fw_q_type: SLURM
rocket_launch: vibes fireworks rlaunch singleshot
nodes: 1
ntasks_per_node: NUMBER OF CORES PER CPU
walltime: MAXIMUM WALL TIME
queue: BATCHING PARTITION
account: ACCOUNT TO TAKE CPU TIME FROM
job_name: DEFAULT JOB NAME
logdir: LOG_FILES
pre_rocket: COMMAND TO RUN BEFORE EACH JOB (if None use null)
post_rocket: COMMAND TO RUN AFTER EACH JOB (if None use null)
```
Finally you can also specify a FWorker by creating an my_fworker.yaml file in \$FW_CONFIG, for example:
```
name: my first fireworker
category: ''
query: '{}'
```
Once all files are in your \$FW_CONFIG folder edit PYTHON_SITE_PACKAGES_DIRECTORY/fireworks/fw_config.py to reflect where to find the correct FireWorks configuration files. The relevant portion changes to the file  should look like this with $FW_CONFIG replaced with the correct path
```
LAUNCHPAD_LOC = $FW_CONFIG/my_launchpad.yaml  # where to find the my_launchpad.yaml file
FWORKER_LOC = $FW_CONFIG/my_fworker.yaml  # where to find the my_fworker.yaml file
QUEUEADAPTER_LOC = $FW_CONFIG/my_qadapter.yaml  # where to find the my_qadapter.yaml file

CONFIG_FILE_DIR = $FW_CONFIG
```
To find where the site-packages file is located run `python -m site` and it should appear in the returned list.


CC: This VIBES specific section should be a separate section in HT, before Phonon Calculations, and should be
called Configuring FHI-vibes for interaction with fireworks

## Creating a `.fireworksrc` file
Now that FireWorks is installed and setup correctly vibes needs to be configured correctly in order to use it effectively.
This file sets up default values for running FireWorks utilities that will be consistent throughout all calculation.
To make this file run 
`vibes template fireworks_configuration` 
CC: The correct command seems to be  vibes template configuration fireworks > ~/.fireworksrc 
CC: Please always put commands in a single line for easier copy/paste
which will create a template `fireworksrc` in the current working directory with the following input:
```
[fireworks]
config_dir: "~/.fireworks"
tasks2queue = ["vibes.relaxation.bfgs.relax", "vibes.fireworks.tasks.calculate_wrapper.wrap_calc_socket", "vibes.k_grid.converge_kgrid.converge_kgrid", "vibes.fireworks.tasks.calculate_wrapper.wrap_calculate", "vibes.fireworks.tasks.md.run" ]

[fireworks.remote]
host = ["remote.host.path"]
config_dir = ["/path/to/remote/home/.fireworks/"]
launch_dir = "."

[fireworks.remote.authorization]
user = remote_username
password = null

[fireworks.remote.launch]
njobs_queue = 0
njobs_block = 500
reserve = True
nlaunches = 0
sleep_time = 60
```
For a complete description of each of these parameters see the [full documentation](../../Documentation/1_general_high_throughput).
Once completed copy `fireworksrc` to `$HOME/.fireworksrc`.
CC: Since not all parameters need to be checked/changed, I would explain the parameters here, even as inline comments,e.g.,
[fireworks]
config_dir: "~/.fireworks"   # Directory containing the *yaml filed
tasks2queue = ["vibes.relaxation.bfgs.relax", "vibes.fireworks.tasks.calculate_wrapper.wrap_calc_socket", "vibes.k_grid.converge_kgrid.converge_kgrid", "vibes.fireworks.tasks.calculate_wrapper.wrap_calculate", "vibes.fireworks.tasks.md.run" ]  ## DO NOT CHANGE
etc...

## Testing if it works
Now that your FireWorks installation should be working properly test it with the vibes FireWorks test in `test/fireworks/test_fireworks.py`.
If the test runs successfully then when you run `lpad get_wflows` you should get the following output (created on should match today's date/time):
```
{
    "state": "COMPLETED",
    "name": "Ni_6d2a2be5a5c1c4549639b55c5403b438b3b0ccf7--1",
    "created_on": "2020-03-13T11:55:30.357000",
    "states_list": "C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C-C"
}

```
If you see that you have successfully set up the high-throughput portions of vibes. To use this on clusters you need to repeat the steps in Installing/Testing FireWorks on each machine you plan to use it on.
