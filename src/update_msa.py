from pathlib import Path


def update_msa(output_path: Path = None):
    """Hack to download latest ms_analysis.py from github GunnerLab/Stable-MCCE repo
    as no updated conda package on channel 'newbooks' exists as of 12-10-2023.
    Note: To update `ms_analysis.py`: back it up or delete it, then rerun: `python update_msa.py`.
    """

    try:
        import ms_analysis as msa

    except ModuleNotFoundError:
        if output_path is None:
            output_path = "./ms_analysis.py"
        else:
            outputdir = str(output_path.joinpath("ms_analysis.py"))
        latest_msa_url = "https://raw.githubusercontent.com/GunnerLab/Stable-MCCE/master/bin/ms_analysis.py"
        curl_cmd = f"curl {latest_msa_url} -o {output_path}".split()

        import subprocess

        try:
            subprocess.run(curl_cmd, check=True)
        except subprocess.CalledProcessError as E:
            raise EnvironmentError(
                f"Could not download ms_analysis.py from 'GunnerLab/Stable-MCCE'.\n{E}"
            )


if __name__ == "__main__":
    update_msa()
