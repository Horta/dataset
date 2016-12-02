f = h5py.File('../../../1000G/hdf5/1000G_stage1.hdf5', 'r')
chrom = 17

sample_ids0 = f[u'geauvadis_variants/chr%s/row_header/sample_ID' % chrom][:]
sample_ids0 = filter(lambda x: x != '0', sample_ids0.tolist())

df0 = df.loc[('ENSG00000000419', 2)]

sample_ids1 = df0['assay'].apply(lambda x: x.split('.')[0]).tolist()

f.close()

sample_ids0_map = dict(zip(sample_ids0, range(len(sample_ids0))))

indices = [sample_ids0_map[s] for s in sample_ids1]
