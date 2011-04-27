import MySQLdb
import MySQLdb.cursors
import csv
import pickle
from itertools import izip_longest
class coge_urls():
    "a class for creating coge links from mysql tables"
    def __init__(self, database, table_name, path):
        "sets up mysql cursor, creates pickle addes to path dir and adds link column to table"
        self.database = database
        db = MySQLdb.connect(host="127.0.0.1", user="root", db = self.database)
        self.cursor=db.cursor(cursorclass=MySQLdb.cursors.DictCursor)
        self.table_name = table_name
        self.cursor=db.cursor(cursorclass=MySQLdb.cursors.DictCursor)
        self.path = path
        
        #check_if_col_exists(table_name, "link2")
        # stmt = "ALTER table {0} add column link varchar(100000)".format(self.table_name) 
        # self.cursor.execute(stmt)
        dbQuery = "SELECT * FROM {0}".format(self.table_name)
        self.cursor.execute(dbQuery)
        table = self.cursor.fetchall()
        file_name = open("{2}/{0}_{1}.pck".format(self.database, self.table_name, self.path) , "w")
        pickle.dump(table, file_name)
        file_name.close()
        
    #def check_if_col_exists():
        
    def grouper(self, n, iterable, fillvalue=None):
        "grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx"
        args = [iter(iterable)] * n
        return izip_longest(fillvalue=fillvalue, *args)
        
    def postion_url(self, disg_ids, padding, blast='n', *col_list):
        "lines up coge based on the cns postion"
        base = "http://synteny.cnr.berkeley.edu/CoGe/GEvo.pl?prog=blast{0}&autogo=1&".format(blast)
        end = 'num_seqs={0};hsp_overlap_limit=0;hsp_size_limit=0'.format((len(col_list))/2)
        pck_file = open('{2}/{0}_{1}.pck'.format(self.database, self.table_name, self.path), "r")
        table_dict = pickle.load(pck_file)
        for row in table_dict:
            row_url = []
            col_tuple = [x for x in self.grouper(2,col_list)]
            for c, p in col_tuple:
                chrm = row[c]
                pos = row[p]
                pos_number = (col_tuple.index((c,p)) + 1)
                dsgid = disg_ids[col_tuple.index((c,p))]
                pos_url = 'dsgid{0}={1}&chr{0}={2}&x{0}={3}&dr{0}up={4}&dr{0}down={4}&'.format(pos_number, dsgid, chrm, pos, padding)
                row_url.append(pos_url)
            inside = ''.join(row_url)
            url = base + inside + pos_url + end
            self.import_url_to_mysql(url, row, col_list)

        
    def import_url_to_mysql(self, url, row,  col_list):
        base_stmt = "UPDATE {0} SET link = '{1}' WHERE".format(self.table_name, url)
        where_stmts = []
        for col in col_list:
            format_number = col_list.index(col)
            where_stmt = " {0} = '{1}'".format(col, row[col])
            where_stmts.append(where_stmt)
        inside_stmt = ' AND '.join(where_stmts)
        
        stmt = base_stmt + inside_stmt
        self.cursor.execute(stmt)
    
    def gene_url(self, padding, *col_list):
        base = "http://synteny.cnr.berkeley.edu/CoGe/GEvo.pl?prog=blastn&autogo=1&"
        end = "num_seqs={0};hsp_overlap_limit=0;hsp_size_limit=0".format(len(col_list))
        pck_file = open('{2}/{0}_{1}.pck'.format(self.database, self.table_name, self.path), "r")
        table_dict = pickle.load(pck_file)
        for row in table_dict:
            row_url = [] 
            for col in col_list:
                gene_number = (col_list.index(col) + 1)
                gene_name = row[col]
                gene_url = 'accn{0}={1}&dr{0}up={2}&dr{0}down={2}&'.format(gene_number, gene_name, padding)
                row_url.append(gene_url)
            inside = ''.join(row_url)
            url = base + inside + end
            self.import_url_to_mysql(url, row, col_list)
            
    def gene_pos_url(self, disg_ids, padding, gene, *col_list):
        "lines up coge based on the cns postion"
        base = "http://synteny.cnr.berkeley.edu/CoGe/GEvo.pl?prog=blastn&autogo=1&"
        end = 'num_seqs=3;hsp_overlap_limit=0;hsp_size_limit=0'
        pck_file = open('{2}/{0}_{1}.pck'.format(self.database, self.table_name, self.path), "r")
        table_dict = pickle.load(pck_file)
        for row in table_dict:
            row_url = []
            col_tuple = [x for x in self.grouper(2,col_list)]
            for c, p in col_tuple:
                chrm = row[c]
                pos = row[p]
                pos_number = (col_tuple.index((c,p)) + 1)
                dsgid = disg_ids[col_tuple.index((c,p))]
                pos_url = 'dsgid{0}={1}&chr{0}={2}&x{0}={3}&dr{0}up={4}&dr{0}down={4}&'.format(pos_number, dsgid, chrm, pos, padding)
                row_url.append(pos_url)
            inside = ''.join(row_url)
            gene_url = 'accn3={0}&dr3up={1}&dr3down={1}&'.format(row[gene], padding)
            url = base + inside + gene_url + end
            print url
            self.import_url_to_mysql(url, row, col_list)
    
gina = coge_urls('rice_v6_maize_v2','non_match_genes','/Users/gturco')
#gina.gene_url(28000, 'qaccn', 'saccn', 'qfeat')
#tommy = coge_urls('trobble_cns', 'gina', '/Users/gturco')
#gina.gene_pos_url([95, 9109], 30000, 'm_accn', 'schr', 's_start', 'qchr', 'q_start')
gina.gene_url(26000, 'maize', 'qaccn', 'sorg' )