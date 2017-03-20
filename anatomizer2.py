""" Class AgentAnatomy. 

Provide the structural features of a protein based on information from 
biological knowledge databases.
"""

import os
import re
import requests
import xml.etree.ElementTree as ET
import xml.dom.minidom
import json
from collections import OrderedDict

class AgentAnatomy(object):
    """ 
    Gather structural features about a protein with given
    HGNC Gene Symbol or UniProt Accession Number
    """

    workdir = 'anatomyfiles'

    species = 'homo_sapiens'

    ensemblserv = 'http://rest.ensembl.org'

    def __init__(self, query):
        """
        Look if query matches a unique Ensembl gene.
        If so, initialize an AngentAnatomy instance. Otherwise, abort.
        """

        self.query = query
        os.makedirs(self.workdir, exist_ok=True)
        
        ensemblext = '/xrefs/symbol/%s/%s?' % (self.species, self.query)
        decoded = self._fetch_ensembl(ensemblext)
        genes = []
        for entry in decoded:
            ensid = entry['id']
            if ensid[0:4] == 'ENSG':
                genes.append(ensid)
        if len(genes) == 1:
            self.ensemblgene = genes[0]
            print('Creating instance of AgentAnatomy with Ensembl Gene ID %s.'
                  % self.ensemblgene)
        else:
            print('Could not find unique Ensembl Gene ID. Aborting.')
            exit()


    def _fetch_ensembl(self,ext):
        r = requests.get(self.ensemblserv+ext,
                         headers={ "Content-Type" : "application/json"})
        if not r.ok:
            r.raise_for_status()
            sys.exit()
        return r.json()


    def _get_hgncsymbol(self):
        ensemblext = '/xrefs/id/%s?' % self.ensemblgene
        xreflist = self._fetch_ensembl(ensemblext)
        for xref in xreflist:
            if xref['db_display_name'] == 'HGNC Symbol':
                self.hgncsymbol = xref['display_id']


    def _get_strand(self):
        ensemblext = '/lookup/id/%s?' % self.ensemblgene
        lookgene = self._fetch_ensembl(ensemblext)
        self.strand = lookgene['strand']


    def _get_transcripts(self):
        ensemblext = '/overlap/id/%s?feature=cds' % self.ensemblgene
        cdslist = self._fetch_ensembl(ensemblext)
        deflist = []
        for cds in cdslist:
            if cds['strand'] == self.strand:
                deflist.append( OrderedDict([ 
                                ('Ensembl_transcr', cds['Parent']), 
                                ('Ensembl_protein', cds['protein_id']) ]) )
        # Remove duplicates (set does not work with dictionaries)
        self.ptnlist = []
        for item in deflist:
            if item not in self.ptnlist:
                self.ptnlist.append(item)


    def _get_hgnctranscr(self):
        for i in range(len(self.ptnlist)):
            enst = self.ptnlist[i]['Ensembl_transcr']
            ensemblext = '/xrefs/id/%s?' % enst
            transcrxreflist = self._fetch_ensembl(ensemblext)
            for txref in transcrxreflist:
                if txref['dbname'] == 'HGNC_trans_name':
                    self.ptnlist[i]['Transcript_name'] = txref['primary_id']


    def _get_uniprotids(self):
        for i in range(len(self.ptnlist)):
            ensp = self.ptnlist[i]['Ensembl_protein']
            ensemblext = '/xrefs/id/%s?' % ensp
            protxreflist = self._fetch_ensembl(ensemblext)
            #self.print_json(protxreflist)
            for pxref in protxreflist:
                if pxref['db_display_name'][:9] == 'UniProtKB':
                    self.ptnlist[i]['UniProt_accession'] = pxref['primary_id']
                    ## Optionally show if from Swiss-prot or TrEMBL
                    #self.ptnlist[i]['UniProt_db'] = pxref['db_display_name'][10:]


    def _fetch_uniprotxml(self,uniprotac):
        """ Retrieve UniProt entry from the web in xml format. """
        if ('uniprot%s.xml' % uniprotac) in os.listdir(self.workdir):
            xmlfile = open('%s/uniprot%s.xml' 
                           % (self.workdir, uniprotac),'r')
            uniprot = xmlfile.read()
            print('Using UniProt entry from file %s/uniprot%s.xml.\n' 
                  % (self.workdir, uniprotac))
        else:
            r = requests.get('http://www.uniprot.org/uniprot/%s.xml' 
                             % uniprotac)
            xmlparse = xml.dom.minidom.parseString(r.text) 
            uniprot = xmlparse.toprettyxml(indent="   ",newl='')
            # Write xml to file to avoid download on future uses
            savefile = open('%s/uniprot%s.xml'
                            % (self.workdir, uniprotac),'w')
            savefile.write(uniprot)
            print('Fetched file from http://www.uniprot.org/uniprot/%s.xml.\n'
                  % uniprotac)
        # Removing default namespace to simplify parsing.
        xmlnonamespace = re.sub(r'\sxmlns="[^"]+"', '', uniprot, count=1)
        root = ET.fromstring(xmlnonamespace)
        return root


    def _get_uniprotdupl(self):
        seen = []
        duplicates = set()
        for ptn in self.ptnlist:
            ac = ptn['UniProt_accession']
            if ac in seen:
                duplicates.add(ac)
            else:
                seen.append(ac)
        # Check UniProt to distinguish ENSPs that have a same UniProt AC.
        for unip in list(duplicates):
            uniprotxml = self._fetch_uniprotxml(unip)
            # Check all the ENSTs from ptnlist that have AC "unip".
            for i in range(len(self.ptnlist)):
                if self.ptnlist[i]['UniProt_accession'] == unip:
                    enst = self.ptnlist[i]['Ensembl_transcr']
                    molecule = uniprotxml.find(".//dbReference[@id='%s']/molecule" % enst)
                    self.ptnlist[i]['UniProt_accession'] = molecule.get('id')


    def _get_length(self):
        for i in range(len(self.ptnlist)):
            ensp = self.ptnlist[i]['Ensembl_protein']
            ensemblext = '/lookup/id/%s?' % ensp
            lookptn = self._fetch_ensembl(ensemblext)
            self.ptnlist[i]['length'] = lookptn['length']


    def _get_canon(self):
        """ Get canonical (primary) transcript from APPRIS """
        r = requests.get('http://apprisws.bioinfo.cnio.es:80/rest/exporter/'
                         'id/%s/%s?methods=appris&format=json' 
                          % (self.species, self.ensemblgene) )
        appris = r.json()
        canonlist = []
        for isoform in appris:
            try:
                if isoform['annotation'] == 'Principal Isoform':
                    canonlist.append(isoform['transcript_id'])
            except:
                pass
        canonset = set(canonlist)
        if len(canonset) == 1:
            self.canon = canonlist[0]
        else:
            print('Cannot find unique canonical (primary) transcript')
            exit()
        for i in range(len(self.ptnlist)):
            if self.ptnlist[i]['Ensembl_transcr'] == self.canon:
                self.ptnlist[i]['Canon'] = 'Yes'



    def print_json(self,data):
        print(json.dumps(data, indent=4))
