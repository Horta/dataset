import numpy as np
import h5py
import blosc
import msgpack
import cPickle as pkl
import tempfile
import pandas as pd
import shutil
import sys
from subprocess import call

pd.set_option('display.width', 1000)

def load_introns():
    root = '/hps/nobackup/stegle/users/horta/dataset/intron/quant_splicing/transcript-qtls'
    folder = tempfile.mkdtemp(dir='/dev/shm')
    shutil.copy(root + '/intron_events_filter2_chrom_pos.pkl.blp', folder + '/introns.pkl.blp')
    call('blpk d ' + folder + '/introns.pkl.blp', shell=True)
    df = pd.read_pickle(folder + '/introns.pkl')
    shutil.rmtree(folder)
    return df

def get_gene_map_intron():
    root = '/hps/nobackup/stegle/users/horta/dataset/intron/quant_splicing/transcript-qtls'
    with open(root + '/gene_map_intron_filter2.msg', 'rb') as f:
        return dict(msgpack.unpackb(f.read()))

def get_window(df, h5, gene, intron):

    df0 = df.loc[(gene, intron)]
    chrom = df0['chrom'][0]
    
    pos = h5[u'geauvadis_variants/chr%s/col_header/pos' % chrom][:]

    nbases = int(10**6 / 2)
    
    pos_start, pos_end = df0['pos_start'][0] - nbases, df0['pos_end'][0] + nbases
    
    index_start, index_end = np.searchsorted(pos, pos_start), np.searchsorted(pos, pos_end)
    
    SNPs = h5[u'geauvadis_variants/chr%s/matrix' % chrom][:, index_start:index_end+1]

    sample_ids0 = h5[u'geauvadis_variants/chr%s/row_header/sample_ID' % chrom][:]
    sample_ids0 = filter(lambda x: x != '0', sample_ids0.tolist())

    df0 = df.loc[(gene0, intron0)]

    sample_ids1 = df0['assay'].apply(lambda x: x.split('.')[0]).tolist()

    sample_ids0_map = dict(zip(sample_ids0, range(len(sample_ids0))))

    indices = [sample_ids0_map[s] for s in sample_ids1]
    
    return SNPs[indices,:]

with h5py.File('../../../1000G/hdf5/1000G_stage1.hdf5', 'r') as h5:

    df = load_introns()
    gene_intron = get_gene_map_intron()
    
    gene0 = gene_intron.keys()[0]
    intron0 = gene_intron[gene0][0]
    
    print(get_window(df, h5, gene0, intron0))

