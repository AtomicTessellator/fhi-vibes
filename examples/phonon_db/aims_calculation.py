import subprocess as sp

def run_aims(folder, aims_out, aims_cmd, force=False, verbose=True):
    """ Run an aims calculation if it has not been performed before."""
    try:
        if 'Have a nice day' in (folder / aims_out).read_text() and not force:
            print(f'aims calculation in ' + str(folder) + ' already finished')
        else:
            raise
    except:
        if verbose: print(f'Run aims calculation in ' + str(folder))
        with open(folder / aims_out, 'wb') as f:
            log = sp.run(aims_cmd.split(), cwd=folder, stdout=f)
