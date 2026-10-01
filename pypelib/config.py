from os import environ


def finish_config_setup(config, base_dir):
    """
    Finish setting up the config

      * load an NCBI API key from file, if needed
      * pass some settings to the environment

    config:
        The snakemake config instance, a dict.
    base_dir:
        pathlib.Path to the directory containing the defaults config file.
    """
    if config.get('ncbi_api_key_file'):
        if config.get('ncbi_api_key'):
            print('[WARNING] Ignoring ncbi_api_key_file setting since '
                  'ncbi_api_key is set already')
        else:
            with open(config['ncbi_api_key_file']) as ifile:
                key_txt = ifile.read().strip()
                try:
                    # testing for expected single line of hex code
                    bytes.fromhex(key_txt)
                except ValueError as e:
                    raise RuntimeError(
                        f'[ERROR] expecting a hex code as api key in '
                        f'{ifile.name}'
                    ) from e
                config['ncbi_api_key'] = key_txt
                # some tools/rules want this in the environment
                environ.setdefault('NCBI_API_KEY', key_txt)

    if slurm_status_cmd := config.get('slurm_status_cmd'):
        slurm_status_cmd = base_dir / slurm_status_cmd
        if slurm_status_cmd.is_file():
            environ.setdefault('SLURM_STATUS_CMD', str(slurm_status_cmd))
        else:
            print(f'[WARNING] config.slurm_status_cmd set to "{slurm_status_cmd}" '
                  f'which is not a file')
