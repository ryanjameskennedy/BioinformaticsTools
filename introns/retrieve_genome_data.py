#!/usr/bin/env python

import os
import gzip
import json
import shutil
from ftplib import FTP
from functools import partial

def create_dir(out_fpath):
    folder = "/".join(out_fpath.split("/")[:-1])
    if not os.path.exists(folder):
        os.makedirs(folder)

def unzip_file(filepath):
    if filepath.endswith(".gz"):
        new_filepath = filepath.rstrip(".gz")
        with gzip.open(filepath, 'rb') as fin:
            with open(new_filepath, 'wb') as fout:
                shutil.copyfileobj(fin, fout)
    os.remove(filepath)
    return new_filepath

def _parse_list_line(line, files=[], subdirs=[], links=None):
    dst = None
    if line.startswith("d"):
        dst = subdirs
    elif line.startswith("-"):
        dst = files
    elif line.startswith("l"):
        dst = links
    else:
        raise ValueError("unknown line type %r" % line[:1])
    if dst is None:
        return
    parts = line.split(None, 8)
    name = parts[-1]
    dst.append(name)

class FTPHost(object):
    def __init__(self, ftp_obj):
        self.ftp_obj = ftp_obj

    def walk(self, directory):
        subdirs, files = self.listdir(directory)
        yield directory, subdirs, files
        for subdir in subdirs:
            for x in self.walk(os.path.join(directory, subdir)):
                yield x

    def listdir(self, directory):
        directory = directory.rstrip("/")
        kwargs = dict(files=[], subdirs=[])
        cb = partial(_parse_list_line, **kwargs)
        self.ftp_obj.dir(directory, cb)
        return kwargs["subdirs"], kwargs["files"]

    def download(self, remote_file, out_fpath):
        with open(out_fpath, 'wb') as fout:
            self.ftp_obj.retrbinary(f"RETR {remote_file}", fout.write)
    
    @classmethod
    def connect(cls, ftp_host):
        ftp_obj = FTP(ftp_host)
        ftp_obj.login()
        return cls(ftp_obj)

    def quit(self):
        self.ftp_obj.quit()

def traverse(ftp_host, out_dir, *args):
    fpath_dict = {}
    ftp = FTPHost.connect(ftp_host)
    for arg in args:
        for root, _, files in ftp.walk(arg):
            for filename in files:
                remote_file = os.path.join(root, filename)
                if "/dna/" in remote_file and remote_file.endswith("dna.toplevel.fa.gz"):
                    out_fpath = out_dir + remote_file.replace("fasta/","").replace("dna/","")
                    create_dir(out_fpath)
                    ftp.download(remote_file, out_fpath)
                    out_fpath = unzip_file(out_fpath)
                    species_name = remote_file.split("/")[-1].split(".")[0].replace("_"," ")
                    if species_name not in fpath_dict:
                        fpath_dict[species_name] = {"fasta": out_fpath}
                    else:
                        fpath_dict[species_name]["fasta"] = out_fpath
                elif remote_file.endswith("48.gff3.gz"):
                    out_fpath = out_dir + remote_file.replace("gff3/","")
                    create_dir(out_fpath)
                    ftp.download(remote_file, out_fpath)
                    out_fpath = unzip_file(out_fpath)
                    if species_name not in fpath_dict:
                        fpath_dict[species_name] = {"gff3": out_fpath}
                    else:
                        fpath_dict[species_name]["gff3"] = out_fpath
    ftp.quit()
    return fpath_dict

def main():
    fpath_dict = traverse("ftp.ensemblgenomes.org", "../data", "/pub/fungi/release-48/fasta", "/pub/fungi/release-48/gff3")
    with open("../data/intron_filepaths.json", 'w+', encoding='utf-8') as fout:
        json.dump(fpath_dict, fout)

if __name__ == "__main__":
    main()
