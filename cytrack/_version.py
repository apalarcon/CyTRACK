import json
import os

pathpkg = os.path.dirname(__file__)
version_file = pathpkg+"/VERSION"
VERSION = open(version_file).read().strip()

update_file = pathpkg+"/LAST_UPDATE"
LAST_UPDATE = open(update_file).read().strip()




version_json = '''
{
 "last_update": "''' + str(LAST_UPDATE) + '''",
 "dirty": false,
 "error": null,
 "contact": "albenis.perez.alarcon@uvigo.es & jose.carlos.fernandez.alvarez@uvigo.es",
 "version": "''' + str(VERSION) + '''",
 "author":"Albenis Pérez-Alarcón, Patricia Coll-Hidalgo, Raquel Nieto, Luis Gimeno"
}
'''  # END VERSION_JSON


def get_versions():
	return json.loads(version_json)
