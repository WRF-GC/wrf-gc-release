.. _hco-sa-download:

########################
Download the source code
########################

The :program:`HEMCO` source code may be downloaded (aka "cloned") with
Git.  By default the :command:`git clone` command will give you the
**main** branch by default:
default.

.. code-block:: console

   $ git clone https://github.com/geoschem/hemco.git HEMCO
   $ cd HEMCO

If you would like a different version of HEMCO you can check out a
different branch.  For example, to check out the **dev** branch, type:

.. code-block:: console

  $ git checkout dev

You can also check out the :program:`HEMCO` source code at the
position of any tag.  For example, to use the :program:`HEMCO` version
3.0.0 code (which is by now an old version), type:

.. code-block:: console

   $ git checkout tags/3.0.0

If you have any unsaved changes, make sure you commit those to a
branch prior to updating versions.
