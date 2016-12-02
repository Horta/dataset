def get_gene_map_intron():
    import msgpack
    import blosc
    root = '/hps/nobackup/stegle/users/horta/dataset/intron/quant_splicing/transcript-qtls'
    with open(root + '/gene_map_intron_filter2.msg', 'rb') as f:
        return dict(msgpack.unpackb(f.read()))
