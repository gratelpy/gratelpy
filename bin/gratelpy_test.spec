# -*- mode: python -*-
a = Analysis(['gratelpy_test'],
             pathex=['.'],
             hiddenimports=[],
             hookspath=None,
             runtime_hooks=None)
a.datas += [('test_reversible_substrate_inhibition.dat', '../gratelpy/tests/test_reversible_substrate_inhibition.dat', 'DATA'),
            ('reversible_substrate_inhibition.txt', '../gratelpy/mechanisms/reversible_substrate_inhibition.txt', 'DATAT')]

pyz = PYZ(a.pure)
exe = EXE(pyz,
          a.scripts,
          exclude_binaries=True,
          name='gratelpy_test',
          debug=False,
          strip=None,
          upx=True,
          console=True )
coll = COLLECT(exe,
               a.binaries,
               a.zipfiles,
               a.datas,
               strip=None,
               upx=True,
               name='gratelpy_test')
