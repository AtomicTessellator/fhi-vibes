**Setup of FireWorks on Computational Resources**
* Connection to the MongoDB virtual machine (VM)
  * Location: tpfwdb.esc.rzg.mpg.de (130.183.206.193)
  * Administrator: Thomas Purcell
  * Accessible via a subset of machines
  	* gatezero.rzg.mpg.de
  	* con01.rzg.mpg.de
  	* Need to have an ssh-key from those servers to get access
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
		  	* openafs-client
		  	* openafs-krb5
		  	* libpam-krb5
		  	* libpam-afs-session
        * once krb5-kds is installed and the /etc/krb5.conf file is updated run kdb5_util create -s
          then run sudo dpkg-reconfigure krb5-kdc
      * Add: ipp-garching.mpg.de mpcdf.mpg.de to CellAlias file
  		* PAM set up is optional
  * Initializing token
  	* kinit username@IPP-GARCHING.MPG.DE may ask for a password here, if no key-chain is used
  	* connection should now be passwordless

* Connecting to the clusters via FireWorks
  * Need fabric 2 and paramiko for remote connections to work
  * Need python-gssapi package to connect with Kerberos (NOT gssapi package)
    * If you are not using a Linux Machine this maybe different
  * Giving hostnames to FireWorks
  	* Always give full name so not draco.mpcdf.mpg.de but draco01.mpcdf.mpg.de

* PostgreSQL Server:
  * VM Machine credentials
    * username: postgres
    * password: HiLDe_DB
  * DB Credentials:
    * user: hilde
    * password: hilde
    * database: phonopy_db

