**Setup of FireWorks on Computational Resources**
* Connection to the MongoDB virtual machine (VM)
  * Location: tpfwdb.esc.rzg.mpg.de (130.183.206.193)
  * Administrator: Thomas Purcell
  * Accessible via a subset of machines
  	* gatezero.rzg.mpg.de
  	* con01.rzg.mpg.de
  	* Need to have an ssh-key from those servers to get access
  	* For an account speak to the administrator
  * MongoDB port: 27017
  	* Accessible from the mpcdf clusters only
  	* To access from your local machine use an ssh tunnel
  	  * ssh -f -N -L local_port:VM_Host:port_on_VM username@mpcdf_cluster_host
  	  * Example: ssh -f -N -L 27017:tpfwdb.esc.rzg.mpg.de:27017 tpurcell@draco.mpcdf.mpg.de

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
  		* Add: ipp-garching.mpg.de mpcdf.mpg.de to CellAlias file
  		* PAM set up is optional
  * Initializing token
  	* kinit username@IPP-GARCHING.MPG.DE may ask for a password here, if no keychain is used
  	* connection should now be passwordless
  * Giving hostnames to FireWorks
  	* Always give full name so not draco.mpcdf.mpg.de but draco01.mpcdf.mpg.de

