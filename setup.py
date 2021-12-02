#!/bin/env/python
#! -*- coding: utf-8 -*-
import setuptools

with open("README.md","r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="SBB",
    version="0.0.19",
    author="Simon Bolduc Beaudoin",
    author_email="Simon.Bolduc.Beaudoin@usherbrooke.ca",
    long_description=long_description,
    long_description_content_type="text/markdown",
    classifiers=[
        "Programming Language :: Python :: 2-3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    packages=setuptools.find_packages(),
    include_package_data=True,
    zip_safe=False, 
    version_config={
        "template": "{tag}",
        "dev_template": "{tag}.post{ccount}+git.{sha}",
        "dirty_template": "{tag}.post{ccount}+git.{sha}.dirty",
        "starting_version": "0.0.19",
        "version_callback": None,
        "version_file": None,
        "count_commits_from_version_file": False,
        "branch_formatter": None,
        "sort_by": None,
        }, 
    setup_requires=["setuptools-git-versioning"],
)
