def load_introns():
    import blosc
    import cPickle as pkl
    import tempfile
    import pandas as pd
    import shutil
    import sys
    from subprocess import call
    root = '/hps/nobackup/stegle/users/horta/dataset/intron/quant_splicing/transcript-qtls'
    
    folder = tempfile.mkdtemp(dir='/dev/shm')
    shutil.copy(root + '/intron_events_filter2_chrom_pos.pkl.blp', folder + '/introns.pkl.blp')
    call('blpk d ' + folder + '/introns.pkl.blp', shell=True)
    df = pd.read_pickle(folder + '/introns.pkl')
    shutil.rmtree(folder)
    return df
