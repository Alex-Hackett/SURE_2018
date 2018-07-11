#!/bin/env python

''' 
 Author: Alexander Hackett 
 Student Number  15323791 
 ahackett@tcd.ie 
 Created for Prof. Archers assignment, Ising Model Simulation 
'''
'''
This is a module that allows sending an email with an arbitrary attachment
It will be used to send the output file from the ising model simulation
to myself to allow analysis on another computer
If run as a programme, will send a test email to myself
Created by Alexander Hackett
15323791


'''
'''
Code largely derivative of the methodology presented in this answer
https://stackoverflow.com/questions/6270782/how-to-send-an-email-with-python
'''
import os
import email
import email.encoders
import email.mime.text
import smtplib
from email.mime import multipart



def email_me(my_email, my_passw, recipients, subject, message, file_name):
    # build the message
    msg = email.mime.multipart.MIMEMultipart()
    msg['From'] = my_email
    msg['To'] = ', '.join(recipients)
    msg['Date'] = email.utils.formatdate(localtime=True)
    msg['Subject'] = subject
    msg.attach(email.mime.text.MIMEText(message))
    
    # build the attachment
    att = email.mime.base.MIMEBase('application', 'octet-stream')
    att.set_payload(open(file_name, 'rb').read())
    email.encoders.encode_base64(att)
    att.add_header('Content-Disposition', 'attachment; filename="%s"' % os.path.basename(file_name))
    msg.attach(att)
    
    # send the message
    srv = smtplib.SMTP('smtp.gmail.com', 587)
    srv.ehlo()
    srv.starttls()
    srv.login(my_email, my_passw)
    srv.sendmail(my_email, recipients, msg.as_string())
    
def main():
    email_me(my_email, my_passw, recipients, subject, message, file_name)
    
if __name__ == '__main__': #If this is run as a test
    #My own details
    my_email = 'alexjanhackett@gmail.com'
    my_passw = '' #Left out for security
    recipients = ['alexjanhackett@gmail.com']
    subject = 'name == main run'
    message = 'This message is being sent because you ran emailer.py as a programme'
    file_name = 'test.txt'#Some test file
    main()

