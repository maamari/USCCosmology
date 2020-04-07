import numpy as np
import matplotlib.pyplot as plt


def ELLIPSE(x0, y0, a, b, theta, n = 1000):

    theta = np.deg2rad(theta)
    t = np.linspace(0, 2*np.pi, n)
    XX = x0 + a*np.cos(theta)*np.cos(t) - b*np.sin(theta)*np.sin(t)
    YY = y0 + a*np.sin(theta)*np.cos(t) + b*np.cos(theta)*np.sin(t)
    e = (XX, YY)

    return e



def Fisher(F, x_fiducial, y_fiducial, xlabel, ylabel, size):

    plt.figure(figsize=size)

    sigma_X, sigma_Y = np.sqrt(np.absolute(np.linalg.inv(F)[0,0])), np.sqrt(np.absolute(np.linalg.inv(F)[1, 1]))
    sigma_XY = np.sqrt(np.absolute(np.linalg.inv(F)[0, 1])) 
    #sigma_X, sigma_Y = np.sqrt(np.abs(F[0, 0])), np.sqrt(np.abs(F[1, 1]))
    #sigma_XY = np.sqrt(np.abs(F[0, 1]))
    
    a2 = np.absolute(0.5 * (sigma_X**2 + sigma_Y**2) + np.sqrt(0.25 * (sigma_X**2 - sigma_Y**2)**2 + sigma_XY**2))
    b2 = np.absolute(0.5 * (sigma_X**2 + sigma_Y**2) - np.sqrt(0.25 * (sigma_X**2 - sigma_Y**2)**2 + sigma_XY**2))
    theta = np.rad2deg(0.5 * np.arctan((2*sigma_XY)/(sigma_X**2 - sigma_Y**2)))

    COLORS = ['C0', 'C1', 'royalblue']
    COLORS0 = ['turquoise', 'c', 'teal']
    COLORS1 = ['orangered', 'orange', 'tomato']
    
    for i, Alpha in enumerate([1.52, 2.3]): # 2,3-sigma 2.48, 3.44

       	E = ELLIPSE(x_fiducial, y_fiducial, Alpha*np.sqrt(a2)/2., Alpha*np.sqrt(b2)/2., theta)

        plt.plot(E[0], E[1], label = str(i+1) + r'$\sigma$', lw=3.0, alpha = 1, color = COLORS[i], zorder = i)
        plt.xlabel(xlabel, size = 22)
        plt.ylabel(ylabel, size = 22)
        plt.legend(prop={'size': 14}, loc = 'best')
        
	#plt.xlim([x_fiducial - x_fiducial*0.0001, x_fiducial + x_fiducial*0.0001])
        #plt.ylim([y_fiducial - y_fiducial*0.0001, y_fiducial + y_fiducial*0.0001])


def Fisher_compare(F0, F1, x_fiducial, y_fiducial, xlabel, ylabel, Alpha, size = (9, 6)):

    plt.figure(figsize=size)

    sigma_X0, sigma_Y0 = np.sqrt(np.absolute(np.linalg.inv(F0)[0,0])), np.sqrt(np.absolute(np.linalg.inv(F0)[1, 1]))
    sigma_XY0 = np.sqrt(np.absolute(np.linalg.inv(F0)[0, 1])) 
    
    sigma_X1, sigma_Y1 = np.sqrt(np.absolute(np.linalg.inv(F1)[0,0])), np.sqrt(np.absolute(np.linalg.inv(F1)[1, 1]))
    sigma_XY1 = np.sqrt(np.absolute(np.linalg.inv(F1)[0, 1])) 
    
    a2_0 = np.absolute(0.5 * (sigma_X0**2 + sigma_Y0**2) + np.sqrt(0.25 * (sigma_X0**2 - sigma_Y0**2)**2 + sigma_XY0**2))
    b2_0 = np.absolute(0.5 * (sigma_X0**2 + sigma_Y0**2) - np.sqrt(0.25 * (sigma_X0**2 - sigma_Y0**2)**2 + sigma_XY0**2))
    theta_0 = np.rad2deg(0.5 * np.arctan((2*sigma_XY0)/(sigma_X0**2 - sigma_Y0**2)))
    
    a2_1 = np.absolute(0.5 * (sigma_X1**2 + sigma_Y1**2) + np.sqrt(0.25 * (sigma_X1**2 - sigma_Y1**2)**2 + sigma_XY1**2))
    b2_1 = np.absolute(0.5 * (sigma_X1**2 + sigma_Y1**2) - np.sqrt(0.25 * (sigma_X1**2 - sigma_Y1**2)**2 + sigma_XY1**2))
    theta_1 = np.rad2deg(0.5 * np.arctan((2*sigma_XY1)/(sigma_X1**2 - sigma_Y1**2)))
    
    E0 = ELLIPSE(x_fiducial, y_fiducial, Alpha*np.sqrt(a2_0)/2., Alpha*np.sqrt(b2_0)/2., theta_0)
    E1 = ELLIPSE(x_fiducial, y_fiducial, Alpha*np.sqrt(a2_1)/2., Alpha*np.sqrt(b2_1)/2., theta_1)

    plt.plot(E0[0], E0[1], c = 'r', label = '50 Mpc cut')
    plt.plot(E1[0], E1[1], c = 'b', label = '10 Mpc cut')
    plt.scatter(x_fiducial, y_fiducial, label = 'Fiducial Value', c = 'k', lw = 0.25)
    plt.xlabel(xlabel, size = 22)
    plt.ylabel(ylabel, size = 22)
    plt.legend(prop={'size':13})
    
