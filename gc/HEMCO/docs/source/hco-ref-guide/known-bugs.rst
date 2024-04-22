.. _hco-known-bugs:

#####################
Known bugs and issues
#####################

Please see our `HEMCO issue tracker on Github
<https://github.com/geoschem/HEMCO/issues>`_ for a list of recent
HEMCO bugs and fixes.

===================
Current bug reports
===================

These `bug reports (listed on the  HEMCO issue tracker)
<https://github.com/geoschem/HEMCO/issues?q=is%3Aissue+is%3Aopen+label%3Abug>`_
are currently unresolved.  We hope to fix these in
future HEMCO releases.

=====================================
Masks cannot be applied to extensions
=====================================

It is currently not possible to :ref:`geographically tag emissions
<cfg-ex-mask-tagged>` computed by :ref:`hco-ext` in the same way that
you would do for :ref:`base emissions <hco-cfg-base>`.  We hope to add
this feature into a future HEMCO release.

====================================================
HEMCO may not recognize alternate spellings of units
====================================================

If a unit string (e.g. :literal:`kg/m2/s`) read from a netCDF
file matches the unit string listed under the :option:`SrcUnit` column
of :ref:`the HEMCO configuration file <hco-cfg>`, then no unit
conversion will happen.

But if the unit string in the file is e.g. :literal:`kg m-2 s-1` and
the unit in the configuration file is :literal:`kg/m2/s`, then HEMCO
detects this as a difference in units, and will try to apply an
automatic conversion that is really unnecssary.

Therefore, we recommend not to rely on HEMCO's automatic unit
capability, and to specfiy all scale factors for unit conversions
explicitly in the configuration file.
