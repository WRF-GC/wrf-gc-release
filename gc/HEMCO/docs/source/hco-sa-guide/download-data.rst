.. _hco-sa-download-data:

###################
Download input data
###################

Before starting a HEMCO standalone simulation, make sure that all of
the relevant emissions and meteorology that you will need for your
simulation are present on disk.

If you are located at an institution where there are several other
HEMCO and/or `GEOS-Chem <https://geos-chem.readthedocs.io>`_ users,
then data for HEMCO standalone might already be located in a shared
folder.  Ask your sysadmin or IT staff.

If you are using HEMCO standalone on the Amazon Web Services EC2
cloud computing platform, then you will have access to an S3 bucket
(:file:`s3://gcgrid/`) with emissions inventories and meteorological data.

If you still need to download data for your HEMCO standalone
simulation, we recommend using the :program:`bashdatacatalog` tool.
For more information, please see our Supplemental Guide entitled
:ref:`bashdatacatalog`.
