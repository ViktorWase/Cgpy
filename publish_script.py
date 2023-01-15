import os
import shutil

def find_version():
	file = open("setup.py", "r")
	text = file.read()
	out = text.split("setuptools.setup(")
	text = out[1]
	out = text.split("version=")
	text = out[1]
	out = text.split(",")
	text = out[0]
	text = text.replace('"', '')
	text = text.replace("'", '')
	file.close()
	return text

def update_version(new_vers, old_vers):
	file = open("setup.py", "r")
	text = file.read()
	file.close()
	file = open("setup.py", "w")
	text = text.replace("version="+'"'+old_vers+'"', "version="+'"'+new_vers+'"')
	file.write(text)
	print(text)
	file.close()


os.system("python -m pip install --user --upgrade setuptools wheel")
os.system("python -m pip install --user --upgrade twine")

version = find_version()
print("The current version is", version)
new_vers = input("Enter next version: ")
update_version(new_vers, version)
os.system("python setup.py sdist bdist_wheel")
#os.system("python -m twine upload --repository-url https://test.pypi.org/legacy/ dist/*")
os.system("python -m twine upload dist/*")

os.system("python -m pip install --user --upgrade pdoc")

os.chdir("cgpy")
os.system("pdoc --html cgp")
os.chdir("..")

cwd = os.getcwd()
source = cwd+"\\"+"cgpy"
files = ["cgp.m.html", "index.html"]
#files = ["cgp.m.html", "index.html", "operation.m.html"]

destination = cwd+"\\"+"docs"
for f in files:
	shutil.move(source+"\\"+f, destination+"\\"+f)


print("Don't forget to commit and push the contents of the docs file, otherwise the documentation won't be updated")

