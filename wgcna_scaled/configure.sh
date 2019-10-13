# clear all acquired locks
sudo rm /var/lib/dpkg/lock
sudo rm /var/lib/dpkg/lock-frontend
sudo rm /var/cache/apt/archives/lock

R_exists=$(which R)
echo $R_exists

if [[ -z "$R_exists" ]]
then
	echo "Installing R......"
	sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
	sudo add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/'
	yes | sudo apt update
	yes | sudo apt install r-base
	echo "R successfully installed"
else
	echo "R is already installed"
fi

# install libopenblas
yes | sudo apt-get install libopenblas-base

echo "lib-open-blas installed"

# install libxml2-dev
yes | sudo apt-get install libxml2-dev
echo "libxml2-dev installed"

yes | sudo apt-get install libssl-dev
echo "libssl-dev installed"

# install libcurl4-openssl-dev
yes | sudo apt-get install libcurl4-openssl-dev
echo "libcurl4-openssl-dev installed"


# install R CRAN packages
sudo su - -c "R -e \"install.packages('devtools', repos='http://cran.rstudio.com/')\""
sudo su - -c "R -e \"install.packages('Rcpp', repos='http://cran.rstudio.com/')\""
sudo su - -c "R -e \"install.packages('RcppEigen', repos='http://cran.rstudio.com/')\""
sudo su - -c "R -e \"install.packages('logging', repos='http://cran.rstudio.com/')\""


echo "R packages installed"

# install R Bioconductor packages
sudo su - -c "R -e \"install.packages('BiocManager'); BiocManager::install('BiocParallel')\""
sudo su - -c "R -e \"BiocManager::install('Biobase')\""
sudo su - -c "R -e \"BiocManager::install('GEOquery')\""

echo "Bioconductor packages installed"

# install WGCNA dependencies
sudo su - -c "R -e \"BiocManager::install('preprocessCore')\""
sudo su - -c "R -e \"BiocManager::install('impute')\""
# remove sqlite lock from instance
sudo rm -r /usr/local/lib/R/site-library/00LOCK-RSQLite
sudo su - -c "R -e \"install.packages('RSQLite', repos='http://cran.rstudio.com/', INSTALL_opts = c('--no-lock'))\""
sudo su - -c "R -e \"BiocManager::install('GO.db')\""
sudo su - -c "R -e \"BiocManager::install('AnnotationDbi')\""
sudo su - -c "R -e \"BiocManager::install('org.Hs.eg.db')\""
sudo su - -c "R -e \"BiocManager::install('org.Mm.eg.db')\""

echo "WGCNA dependencies installed"

# install WGCNA package 
sudo su - -c "R -e \"install.packages('WGCNA', repos='http://cran.rstudio.com/')\""

echo "WGCNA installed"

# install AWS CLI
yes | sudo apt install python-pip
sudo -H pip install --upgrade awscli

echo "AWS CLI installed"
