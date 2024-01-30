#!/usr/bin/python3
'''
Download and install libvna in Native Windows
'''

import argparse
import json
import os
import re
import requests
import shutil
import subprocess
import sys
import tempfile
import urllib.request

progname = os.path.basename(sys.argv[0])

URL = 'https://api.github.com/repos/scott-guthridge/libvna/releases/latest'


# libvna-0.3.9-windows-x64-stdc.zip
# libvna-0.3.9-windows-x64-msvc.zip

def main():
    # Parse Options
    argp = argparse.ArgumentParser(description='Download and install libvna')
    argp.add_argument('-A', '--add-path',
                      choices=['no', 'user', 'global'],
                      default='U', help='Control path addition')
    argp.add_argument('-M', '--MSVC', action='store_true',
                      help='Use MSVC complex numbers')
    argp.add_argument('-p', '--pkgconfig-prefix', type=str,
                      help='set path in pkg-config files')
    argp.add_argument('install_directory', type=str, nargs='?',
                      help='install directory path')
    args = argp.parse_args()

    if args.MSVC:
        suffix = '-windows-x64-msvc.zip'
        default_target = 'C:\\Program Files\\libvna-msvc'
    else:
        suffix = '-windows-x64-stdc.zip'
        default_target = 'C:\\Program Files\\libvna-stdc'

    if args.install_directory:
        target_dir = args.install_directory
    else:
        target_dir = default_target

    if args.pkgconfig_prefix:
        pkgconfig_prefix = args.pkgconfig_prefix
    else:
        pkgconfig_prefix = re.sub(r'\\', r'/', target_dir)

    # Find the latest windows zip file.
    response = requests.get(URL)
    if response.status_code != 200:
        print(f'{progname}: error: {URL}: status {response.status_code}')
        exit(1)
    data = response.json()
    download_url = None
    zip_filename = None
    for asset in data['assets']:
        if 'browser_download_url' in asset:
            download_url = asset['browser_download_url']
            if download_url.endswith(suffix):
                zip_filename = os.path.basename(download_url)
                break
    else:
        print(f'{progname}: libvna-*-{suffix}: not found in assets',
              file=sys.stderr)
        exit(3)

    # Download and extract the zip file
    temp_dir = tempfile.mkdtemp()
    try:
        zip_path = os.path.join(temp_dir, zip_filename)
        urllib.request.urlretrieve(download_url, zip_path)
        os.makedirs(target_dir, exist_ok=True)
        subprocess.run(['tar', '-xf', zip_path, '-C', target_dir])
    finally:
        shutil.rmtree(temp_dir)

    # Repair the pkg-config files.
    replacements = [
        ['prefix=', pkgconfig_prefix],
        ['exec_prefix=', pkgconfig_prefix],
        ['lib_dir=', '/'.join([pkgconfig_prefix, 'lib'])],
        ['include_dir=', '/'.join([pkgconfig_prefix, 'include'])],
    ]
    for pc_filename in ['libvna.pc', 'yaml-0.1.pc']:
        pc_pathname = os.path.join(target_dir, 'lib', 'pkgconfig',
                                   pc_filename)
        new_contents = ''
        with open(pc_pathname, 'r+') as pc_file:
            for line in pc_file:
                for repl in replacements:
                    if line.startswith(repl[0]):
                        line = repl[0] + repl[1] + '\n'
                        break
                new_contents += line

            # Replace the file contents.
            pc_file.seek(0)
            pc_file.write(new_contents)

    # Add the bin directory to PATH
    if args.add_path != 'no':
        cmd = ['setx']
        if args.add_path == 'global':
            cmd.append('/M')
        new_path = os.environ['PATH']
        if not new_path.endswith(';'):
            new_path += ';'
        new_path += os.path.join(target_dir, 'bin') + ';'
        cmd += ['PATH', new_path]
        if os.environ['PATH'] and target_dir not in os.environ['PATH']:
            subprocess.run(cmd)


try:
    main()
except Exception as e:
    print(f'{progname}: exception: {e}', file=sys.stderr)
    sys.exit(99)
