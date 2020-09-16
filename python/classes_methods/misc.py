#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 27 13:44:31 2020

Contains functions to run the user menu and the interaction of user and tool.

@author: jonaszbinden, adopted from Matthias Brennwald EAWAG Dübendorf
"""
import sys
import warnings
import time
# try:
# 	from termcolor import colored
# 	do_color_term = True
# except ImportError:
# 	do_color_term = False
# 	print ('*** NOTE: Please install the python termcolor package for colored warning messages! ***')

# check Python version and print warning if we're running version < 3:
if ( sys.version_info[0] < 3 ):
	warnings.warn("ruediPy / misc class is running on Python version < 3. Version 3.0 or newer is recommended!")


class misc:
    print('*** Welcome to the CRIRES+ Observation Planner ***')
    
#     @staticmethod
#     def warnmessage(unit, msg):
#         """
#         c.warnmessage(caller,msg)

# 		Print a warning message
# 		
# 		INPUT:
# 		caller: caller label / name of the calling object (string)
# 		msg: warning message
# 		
# 		OUTPUT:
# 		(none)
#         """
# 		
#         print('\a') # get user attention using the terminal bell
#         M = '***** WARNING from ' + unit + ' at ' + misc.now_string() + ': ' + msg + '\n'
# 		
#         # if do_color_term:
#         # 	print(colored(M,'red'))
#         # else:
#         #     print(M)	


	########################################################################################################
	

    @staticmethod
    def wait_for_enter(msg='Press ENTER to continue.'):
        """
		misc.wait_for_enter(msg='Press ENTER to continue.')
		
		Print a message and wait until the user presses the ENTER key.
		
		INPUT:
		msg (optional): message
		
		OUTPUT:
		(none)
		"""
		
        print ('\a') # get user attention using the terminal bell
		
        input(msg)
        

	########################################################################################################
	

    @staticmethod
    def ask_for_value(msg='Enter value = '):
        """
  		x = misc.ask_for_value(msg='Enter value = ')
  		
  		Print a message asking the user to enter something, wait until the user presses the ENTER key, and return the value.
  		
  		INPUT:
  		msg (optional): message
  		
  		OUTPUT:
  		x: user value (string)
  		"""  		
		
        print ('\a') # get user attention using the terminal bell
		
       	x = input(msg)
        
		
        return x

    @staticmethod
    def user_menu(menu,title='Choose one of the following options'):
        """
		x = misc.user_menu(menu,title='Choose an option')
		
		Show a "menu" for selection of different user options, return user choice based on key pressed by user.
		
		INPUT:
		menu: menu entries (tuple of strings)
		title (optional): title of the menu (default='Choose an option')
		
		OUTPUT:
		x: number of menu choice
		
		EXAMPLE:
		k = misc.user_menu( title='Choose dinner' , menu=('Chicken','Burger','Veggies') )
		"""		
		
        print ('\a') # get user attention using the terminal bell
        N = len(menu)
        do_menu = True
        while do_menu:
			
            print( '\n' + title + ':' )
            for i in range(N):
                print(' ' + str(i+1) + ': ' + menu[i] )
            
            ans = input( 'Enter number: ' )
			
            try:
                ans = int(ans) # try converting from string to integer number
            except ValueError:
                ans = -1
				
            if int(ans) in range(1,N+1):
                do_menu = False

            if do_menu:
                print('\nInvalid input. Try again...')		
		
        return ans							
