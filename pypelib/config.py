from os import environ

from snakemake.common.configfile import load_configfile


DEFAULTS_CONFIG_FILE = 'default.conf'


def finish_config_setup(config):
    """
    Finish setting up the config

      * apply default configuration
      * load an NCBI API key from file, if needed
      * pass some settings to the environment

    config:
        The snakemake config instance, a dict.
    """
    for key, value in load_configfile(DEFAULTS_CONFIG_FILE).items():
        if value is None:
            # treat None as unset
            continue
        if key not in config:
            config[key] = value

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

