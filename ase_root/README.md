# web-page

## Instructions
### Clone the most recent ASE to a separate directory
```bash
$ git clone https://gitlab.com/ase/ase.git
```

### Copy the modified scripts from this repo's `ase_db` folder to the cloned ASE repo:
```bash
$ cp -r /web-page/ase_root/* /ase-3.20.0b1
$ cp -r /web-page/ase_db/* /ase-3.20.0b1/ase/db
```

### Go to cloned ASE repo and (re)install ASE package by
```bash
$ cd ase-3.20.0b1/
$ python setup.py install
```

### Copy all `static` subdirectories to installation:
```bash
cp -r /web-page/ase_db/static/css ~/anaconda3/lib/python3.6/site-packages/ase-3.20.0b1-py3.6.egg/ase/db/static/
cp -r /web-page/ase_db/static/images ~/anaconda3/lib/python3.6/site-packages/ase-3.20.0b1-py3.6.egg/ase/db/static/
cp -r /web-page/ase_db/static/js ~/anaconda3/lib/python3.6/site-packages/ase-3.20.0b1-py3.6.egg/ase/db/static/
```
NOTE: In this current version, these subfolders are not automatically unpacked to the install dir. during the Python setup.

### Copy and unzip the jsmol.zip to the ase/db/static
To see the 3D structure, one needs to copy the jsmol file (`ase_db/static/jsmol.zip`)
to the static folder (`Installation_of_ase/ase/db/static/jsmol`)

An example is given as follows,
```bash
cp -r ase/db/static/jsmol/* ~/anaconda3/lib/python3.6/site-packages/ase-3.20.0b1-py3.6.egg/ase/db/static/jsmol
```
### View database

```bash
$ ase db tp02.db -w
```
Then open browser with the url of http://0.0.0.0:5000/

### Running on CMS
```bash
$ conda activate tpdb
(tpdb) $ cd /scratch/tpdb/github/web-page
(tpdb) $ git pull origin master
(tpdb) $ cp -r ase_root/* /scratch/tpdb/ase-3.20.0b1
(tpdb) $ cp -r ase_db/* /scratch/tpdb/ase-3.20.0b1/ase/db
(tpdb) $ cd /scratch/tpdb/ase-3.20.0b1
(tpdb) $ pip install -e .
(tpdb) $ cp -r /scratch/tpdb/github/web-page/ase_db/static/css ~/miniconda3/envs/tpdb/lib/python3.8/site-packages/ase-3.20.0b1-py3.8.egg/ase/db/static/
(tpdb) $ cp -r /scratch/tpdb/github/web-page/ase_db/static/images ~/miniconda3/envs/tpdb/lib/python3.8/site-packages/ase-3.20.0b1-py3.8.egg/ase/db/static/
(tpdb) $ cp -r /scratch/tpdb/github/web-page/ase_db/static/js ~/miniconda3/envs/tpdb/lib/python3.8/site-packages/ase-3.20.0b1-py3.8.egg/ase/db/static/
(tpdb) $ cd /scratch/tpdb/
(tpdb) $ nohup ase db tp01.db -w &
```
Hit enter once more after the last command.
