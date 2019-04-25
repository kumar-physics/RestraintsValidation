import sys
import wget
import urllib




def get_data_from_gird(flist):
    f = open(flist,'r').read()
    pdb_ids = []
    for r in f.split("\n"):
        pdb_ids+=r.split(" ")
    pl = [i.lower() for i in pdb_ids]
    f = open('idlist.txt','w')
    f.write("\n".join(pl))
    f.close()

    # for pdb in pdb_ids:
    #     try:
    #         url = 'http://www.bmrb.wisc.edu/ftp/pub/bmrb/nmr_pdb_integrated_data/coordinates_restraints_chemshifts' \
    #                 '/all/nmr-star/{}/{}_linked.str'.format(pdb.lower(),pdb.lower())
    #         wget.download(url)
    #         pdb_url = 'https://files.rcsb.org/download/{}.cif'.format(pdb)
    #         wget.download(pdb_url,'{}.cif'.format(pdb.lower()))
    #
    #     except urllib.error.HTTPError:
    #         print (pdb_url)

if __name__ == "__main__":
    get_data_from_gird('list_from_pdb.txt')