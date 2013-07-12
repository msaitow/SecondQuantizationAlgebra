#    file:  sqaOptions.py
#  author:  Eric Neuscamman
#    date:  March 30, 2009
# summary:  Contains options that the user may modify.
#
# (c) 2008-2009 Eric Neuscamman (eric.neuscamman@gmail.com)
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# In addition, any modification or use of this software should
# cite the following paper:
#
#   E. Neuscamman, T. Yanai, and G. K.-L. Chan.
#   J. Chem. Phys. 130, 124102 (2009)


class sqaOptions:
  "A class to hold options for the second quantization algebra program."

  def __init__(self):
    # Set default options
    self.verbose = True
    self.alpha_type = ('alpha',)
    self.beta_type = ('beta',)
    self.core_type = ('core',)
    self.active_type = ('active',)
    self.virtual_type = ('virtual',)

# Create an object of the options class
options = sqaOptions()
