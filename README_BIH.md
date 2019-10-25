
# Setting up your NGS analysis environment for the ´ngs´ repo

## Run ´setup_ngs_env.sh´

The scripts was developed to automatize setup as much as possible.
Run it using default parameters as:

´´´
setup_ngs_env.sh
´´´

This will download and install miniconda, install all of the necessary packages that can be downloaded and installed without limiting license restrictions and informs on how to proceed next. Files that cannot be downloaded by the script will be informed about.

There are four optional parameters that can be passed to this script:

´´´
set_ngs_env.sh \
	new_miniconda_installation_folder \
	new_miniconda_environment_name \
	new_ngs_directory_to_save_all_other_files \
	server_name
´´´

´new_miniconda_installation_folder´ is where miniconda will be installed to; default is ´${HOME}/work/miniconda´.

´new_miniconda_environment_name´ is the name of the miniconda enviroment that will be created and to which all of the necessary packages will be installed; default is ´ngs´. This is the environment that will have to be activated to use these packages.

´new_ngs_directory_to_save_all_other_files´ is, as the name suggests, the folder where all other necessary files which are not miniconda-related will be downloaded to. In this folder, which defaults to the current directory ´.´, a ´NGS´ parent folder will be created to harbor all of these files. On the BIH server, I recommend ´"${HOME}/work/´.

All of the newly installed miniconda and other package toools will be added to ´$PATH´ via an entry in ´.bashrc´. If the ´server_name´ is set to ´BIH´, which is also its default value, these changes to ´$PATH´ are only activated when one is *not* logged in to the login nodes. This can help to prevent starting by accident an analysis directly on the login nodes.

Thus, calling ´setup_ngs_env.sh´ without any specified arguments is analoguous to calling
´´´
setup_ngs_env.sh \
	"${HOME}/work/miniconda" \
	'ngs' \
	'.' \
	'BIH'
´´´

## Completing setup
Follow the instructions given by ´setup_ngs_env.sh´ to incorporate any missing files into the ´NGS´ setup directory.

Note: GATK and ANNOVAR both come with limiting licenses. *Please register at the Broad Institute and ESPECIALLY the ANNOVAR website to comply with these licenses!!!*


## Using / Activating the set up environment
To activate the miniconda environmnent, use:
´´´
conda activate new_miniconda_environment_name
´´´

The bpipe pipelines can then be run as usual. All other scripts / programs are run as usual in an interactive sesssion or as you would run them on the server via qsub.
