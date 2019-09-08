smSPK
=====
Install dependencies:
if using virtual env omit `--user` parts
```bash
pip install -r requirements.txt --user
```
For mosek set environment variable `MOSEKLM_LICENSE_FILE` to your license file
install `mosek` package with
```bash
pip install -f https://download.mosek.com/stable/wheel/index.html Mosek --user
```
