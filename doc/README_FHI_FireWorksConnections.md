**Setup of FireWorks on Computational Resources**
* Connection to the MongoDB virtual machine (VM)
  * Location: tpfwdb.esc.rzg.mpg.de (130.183.206.193)
  * Administrator: Thomas Purcell
  * Accessible via a subset of machines
    * gatezero.rzg.mpg.de
    * con01.rzg.mpg.de
    * Need to have an ssh-key from those servers to get access
      * Better to you local host-specific keys and use an ssh-agent to get access
      * See general_unix_tips channel for an illustration on how to do that
    * For an account speak to the administrator
  * MongoDB port: 27017 (individual ports will be given by the Administrator)
    * Accessible from the mpcdf machines only
    * To access from your local machine using the vpn

* Setting up passwordless login to mpcdf clusters on FHI laptops
  * Kerberos client setup
    * setup instructions: https://www.mpcdf.mpg.de/services/desktop-support/configuring-your-pc
    * modifications
      * All packages available from package manager
        * krb5-kdc
        * libkrb-dev
        * openafs-client
        * openafs-krb5
        * libpam-krb5
        * libpam-afs-session
        * once krb5-kdc is installed and the /etc/krb5.conf file is updated run kdb5_util create -s
          then run sudo dpkg-reconfigure krb5-kdc
      * Add: ipp-garching.mpg.de mpcdf.mpg.de to CellAlias file
      * PAM set up is optional
  * Initializing token
    * kinit username@IPP-GARCHING.MPG.DE may ask for a password here, if no key-chain is used
    * connection should now be passwordless

* Setting up my_qadapater.yaml for a SLURM Account
  * \_fw_name: CommonAdapter # Always this unless you modify the Queue Adapter class
  * \_fw_q_type: SLURM
  * rocket_launch: rlaunch_hilde singleshot
  * nodes: 1
  * ntasks_per_node: 32 # Number of CPUs per node
  * cput_per_task: 1 # Number of Threads per CPU
  * walltime: '00:30:00'
  * queue: express
  * account: null # This should be null on mpcdf (used if you have an account with hour tracking on a cluster)
  * job_name: null
  * logdir: /u/MPCDF_USER/.fireworks/logs/ # Must be absolute path
  * pre_rocket: module load XXXX YYYY ZZZZ
  * post_rocket: null

* Connecting to the clusters via FireWorks
  * Need fabric 2 and paramiko for remote connections to work
  * Need python-gssapi package to connect with Kerberos (NOT gssapi package)
    * If you are not using a Linux Machine this maybe different
  * Giving hostnames to FireWorks
    * Always give full name so not draco.mpcdf.mpg.de but draco01.mpcdf.mpg.de
* Testing your FireWorks Installation
  * Use the script in $HILDE_HOME/tests/fireworks/fireworks_install_test.py
  * If you see *******All Tests Successful******* at the end then it is installed correctly
  * If you see *******Resetting LaunchPad******* at the end
    * Then there is a problem with your local connection to the LaunchPad
    * If Authentication error then my_launchpad.yaml has the wrong username/password
    * If it hangs and you get a pymongo.errors.ServerSelectionTimeoutError: 130.183.206.193:port number: timed out then you can't connect to the server
      * Check your host/port number are correct
      * Check you're on the right VPN
    * If you see *******Testing the local version of HiLDe's FireWorks modifications******* at the end
      * Local HiLDe modifications are not working.
      * Check HiLDe installation
      * Post error in slack channel
    * If you see *******Testing the remote connection with FireWorks******* at the end
      * Problem with remote connection
      * Check remote HiLDe installation
      * Post error in the slack channel
* PostgreSQL Server:
  * VM Machine credentials
    * username: postgres
    * password: HiLDe_DB
  * DB Credentials:
    * user: hilde
    * password: hilde
    * database: phonopy_db

