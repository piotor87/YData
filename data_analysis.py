import sqlite3,csv
import numpy as np
import pandas as ps
import os
cwd = os.getcwd() +'/'
sqlitePath = cwd  +  'events.sqlite'
dataPath = cwd + 'events.csv'
import pylab
import matplotlib as mpl
mpl.use('TkAgg')
from matplotlib import pyplot as plt
from verkko.plots import matplotlibHelperFunction as HF

from lifelines import AalenAdditiveFitter 
from lifelines.statistics import logrank_test
from lifelines import KaplanMeierFitter
from lifelines.plotting import add_at_risk_counts
from lifelines import NelsonAalenFitter
naf = NelsonAalenFitter()
sqlitePath = cwd + '/' +  'events.sqlite'



def plot_life(order = True, N = 20):
    '''
    Plots an example of the typical birth/death mechanics for groups B and C
    '''

    
    pylab.ioff()
    finalFigPath = cwd  + 'lifetime.pdf'  
    print finalFigPath
    
    # SET UP FIG
    fig = HF.setFigure()
    gs = mpl.gridspec.GridSpec(2, 1)

    ax1 = fig.add_subplot(gs[0, 0] )
    ax2 = fig.add_subplot(gs[1, 0] )
    axes = [ax1,ax2]
    

    groups = ['b','c'] 
    for i,ax in enumerate(axes):

        groupName = groups[i]

        #import data
        lifetimes,observed = return_unsub_time_arrays(groupName)
        birthtimes = return_birthdays(groupName)

        #shuffle
        idx = np.random.choice(np.arange(len(lifetimes)), N, replace=False)
        lifetimes = np.array(lifetimes[idx])
        event_observed = np.array(observed[idx])
        birthtimes = birthtimes[idx]
        
        if order:
            """order by length of lifetimes; probably not very informative."""
            ix = np.argsort(lifetimes, 0)
            print ix,lifetimes
            lifetimes = lifetimes[ix]
            event_observed = event_observed[ix]
            birthtimes = birthtimes[ix]

        for i in range(N):
            c = "#A60628" if event_observed[i] else "#348ABD"
            ax.hlines(N - 1 - i, birthtimes[i], birthtimes[i] + lifetimes[i], color=c, lw=1)
            m = "|" if not event_observed[i] else 'o'
            ax.scatter((birthtimes[i]) + lifetimes[i], N - 1 - i, color=c, s=10, marker=m)
            ax.scatter((birthtimes[i]),  N - 1 - i, color=c, s=10, marker='o')

        ax.tick_params(    
            axis='y',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            left='off',      # ticks along the bottom edge are off
            right='off',         # ticks along the top edge are off
            labelleft='off') # labels along the bottom edge are off

        ax.set_xlim([0,lifetimes.max()+1])
    ax2.set_xlabel('Time from beginning of analysis')
   
    
    fig.savefig(finalFigPath)
    plt.close()


                
def plot_survival():
    '''
    Plots the KM and NA estimates
    '''
    pylab.ioff()
    finalFigPath = cwd + 'survival.pdf'  
    print finalFigPath


    # SET UP FIG
    fig = HF.setFigure()
    gs = mpl.gridspec.GridSpec(2, 1)


    ax1 = fig.add_subplot(gs[0, 0] )
    ax2 = fig.add_subplot(gs[1, 0] )
    colors = ['red','blue']
    groups = ['b','c']

    ax1.set_ylabel(r"$\hat{S}(t)$", fontsize = 8)
    ax2.set_ylabel(r"$\hat{\Lambda} (t)$", fontsize = 8)
    ax2.set_xlabel('Days from subscription')

    # lists where I store data needed for the logRank test and the printing of at_risk_counts
    logrankData = []
    kmfList = []
    for i,group in enumerate(groups):

        #import fitter
        kmf = KaplanMeierFitter()
        #import data from group
        durations,observed = return_unsub_time_arrays(group)
        # append data 
        logrankData.append(durations)
        logrankData.append(observed)
        kmfList.append(kmf)

        #fit
        kmf.fit(durations, event_observed = observed, label = 'Group ' +str(group),ci_labels = ['lower bound','upper bound'])
        naf.fit(durations, event_observed = observed, label = 'Group ' +str(group))
        #plot
        naf.plot(ax = ax2, color = colors[i])
        kmf.plot(ax = ax1, color = colors[i])

        
    add_at_risk_counts(kmfList[0],kmfList[1], ax=ax1)

    results = logrank_test(logrankData[0],logrankData[2],logrankData[1],logrankData[3], alpha=.99)

    results.print_summary()
    ax1.set_xlabel('')

    ax2.set_xlabel('Days from subscription')
    ax2.legend(loc = 'lower right')
    ax1.set_ylim([0,1])
    fig.savefig(finalFigPath)
                



def profit(T,groupName = 'b',worstCase = 'True'):

    if  groupName =='b':
        return 5*(0.37 + 0.63*np.ceil(T/7.))
    if groupName =='c':                
        if worstCase is False:
           c =  np.average(np.arange(1,8)*0.8)
           return (0.37*c + 0.63*T*0.8)

        else:
            return 0.8*(0.37+0.63*T)


def plot_time(weeks = 11):
    '''
    Plots an example of the typical birth/death mechanics for groups B and C
    '''
    
    tFunc = np.vectorize(profit)
    
    pylab.ioff()
    finalFigPath = cwd  + 'profit.pdf'  
    print finalFigPath
    
    # SET UP FIG
    fig = HF.setFigure()
    gs = mpl.gridspec.GridSpec(1, 1)

    ax = fig.add_subplot(gs[0, 0] )
    
    time = np.arange(weeks*7)+1
    profitB = tFunc(time,'b')
    profitC = tFunc(time,'c',worstCase = True)
    profitBestC= tFunc(time,'c',worstCase = False)
          
    ax.step(time,profitB, 'r',label = r'$P_{B}(T)$')
    ax.plot(time,profitC, 'b--', label = r'$P_{C_{\!^{\_}}}(T) $')
    ax.plot(time,profitBestC, 'b+', label = 'P_{C_{+}}(T) ')
    
    ax.fill_between(time, profitC, profitBestC, alpha = 0.1)

    ax.set_xlabel('Time in Weeks')
    ax.set_ylabel('Average Income per User')

    r = np.arange(1,weeks+1)
    ticks_pos = [-2.5+7*i for i in r]
    print ticks_pos

    ticks_labels = r
    plt.xticks([7*i+1 for i in r], ticks_labels)
    ax.legend(loc = 'lower right')
    fig.savefig(finalFigPath)
    plt.close()


def return_birthdays(groupName = 'b'):
    '''
    Returns subscription event in terms of days since the start of the experiment, where t = 0 is set with the first event recorded.
    '''

    cmd = """
    SELECT ROUND(julianday(s.date) - (SELECT (MIN(julianday(date))) FROM events) """ + '' + """ - 0.5)   
    FROM subscribed AS s     
    WHERE uid IN (SELECT uid FROM events WHERE groupName = '""" + str(groupName) + """')"""

    with sqlite3.connect(sqlitePath) as conn:
        cursor =conn.cursor()
        cursor.execute(cmd,)
        data = cursor.fetchall()

    return np.array([int(elem[0]) for elem in data])

    
def return_unsub_time_arrays(groupName = 'a'):

    '''
    Returns two arrays.
    Durations: number of days (rounded up) since between subscription and unsubscription. In case of user still subscribed, the difference is calculated from the time stamp of the last event
    Observed: Boolean array that tells us if there has been a unsubscription or not.
    '''

    cmd = """

    SELECT ROUND( IFNULL(julianday(u.date), (SELECT (MAX(julianday(date))) FROM events) ) - julianday(s.date) + 0.5),
    ROUND( IFNULL( julianday(u.date)/julianday(u.date),0))
    
    FROM subscribed AS s 
    LEFT JOIN unsubscribed AS u USING (uid)
    WHERE uid IN (SELECT uid FROM events WHERE groupName = '""" + str(groupName) + """')"""

    with sqlite3.connect(sqlitePath) as conn:
        cursor =conn.cursor()
        cursor.execute(cmd,)
        data = cursor.fetchall()

        durations = np.array([int(elem[0]) for elem in data])
        observed = np.array([bool(elem[1]) for elem in data])
    return durations,observed




def return_z_test(p1 = 326/332.,p2 = 330/332., n1= 332.,n2 = 332.):
    '''
    Returns z-test score given probabilities and populations
    '''
    n1 = float(n1)
    n2 = float(n2)
    print p1,p2,n1,n2
    p = (n1*p1+n2*p2)/(n1 + n2)
    print p
    z = (p1-p2)/( np.sqrt( p*(1-p)*(1/n1 + 1/n2))   )

    return np.abs(z)
                


# CONVERSION RATE FOR QUESTION 2
def return_conversion_rates():

    '''
    Returns, for each group, the number of users who opened an account, played a song after creating the account and the conversion rate
    '''

    
    with sqlite3.connect(sqlitePath) as conn:
        cursor =conn.cursor()
        for groupName in ['a','b','c']:
            ac = 333. - return_accounts_created(groupName)
            sp = return_songs_played(groupName)

            print groupName, ac,sp,float(sp)/ac
            

def return_accounts_created(groupName = 'a'):

    with sqlite3.connect(sqlitePath) as conn:
        cursor = conn.cursor()

        cmd = """
    
        SELECT COUNT(DISTINCT(uid))
        FROM events  
        WHERE uid 
        NOT IN (SELECT DISTINCT(uid) FROM account) 
        AND groupName = (?)
    """
        cursor.execute(cmd,(groupName))
        data = cursor.fetchall()[0][0]

    return data


def return_songs_played(groupName = 'a'):
    
    with sqlite3.connect(sqlitePath) as conn:
        cursor = conn.cursor()

        cmd = """
        SELECT COUNT(DISTINCT(uid)) 
        FROM account 
        LEFT JOIN events USING (uid)
        WHERE events.date > account.date
        AND event = 'song_played'
        AND groupName = (?)"""
        cursor.execute(cmd,(groupName))

        data = cursor.fetchall()[0][0]

    return data



def return_answer_one():
    '''

    Returns averages for question 1
    '''


    with sqlite3.connect(sqlitePath) as conn:

        cursor = conn.cursor()
        cmd = """
        SELECT AVG(cnt)
        FROM (        
        SELECT count(*) AS cnt
        FROM subscribed 
        LEFT JOIN events USING (uid)
        WHERE (events.date < subscribed.date) AND event in ('minigame_played','song_played') 
        GROUP BY (uid) )      
        """
        cursor.execute(cmd,)
        data = cursor.fetchall()
        print data

        
        cmd = """
        SELECT AVG(cnt)
        FROM (
        SELECT count(e.uid)  AS cnt
        FROM subscribed AS s
        LEFT JOIN events AS e
        ON (e.uid = s.uid) 
        AND (e.date < s.date) 
        AND e.event in ('minigame_played','song_played') 
        GROUP BY (s.uid) )
        """

        cursor.execute(cmd,)
        data = cursor.fetchall()
        print data


        cmd = """
        SELECT COUNT(*) 
        FROM subscribed 
        LEFT JOIN events
        USING (uid)
        WHERE events.date < subscribed.date
        AND  event in ('minigame_played','song_played')
        """
        cursor.execute(cmd,)
        data = cursor.fetchall()
        entries = data[0][0]
        print entries

        cmd = """

        SELECT COUNT(*)
        FROM subscribed
        """
        cursor.execute(cmd,)
    
        subUsers = cursor.fetchall()[0][0]
        print subUsers

        print float(entries)/subUsers



def create_all():

    create_database()
    create_subscribed_table()
    create_unsubscription_table()
    create_account_creation_table()
    create_first_song_played()
def create_database():
    '''
    Functions that creates a sqlite database from the events.
    '''

    tableName = 'events' 

    with sqlite3.connect(sqlitePath) as conn:

        cursor = conn.cursor()
        conn.execute("""DROP TABLE IF EXISTS """ + tableName )
        conn.execute("""CREATE TABLE IF NOT EXISTS """ + tableName +  """ (event TEXT, uid TEXT, groupName TEXT, date DATE) """)


        # READ CSV DATA
        csvReader = csv.reader(open(dataPath), delimiter=',', quotechar='"')
        next(csvReader, None)  # skip the first line

        
        # INSERT DATA IN TABLE
        for row in csvReader:
            conn.execute('INSERT INTO events (event, uid, groupName, date) values (?, ?, ?,?)', row)
            

        # CREATE INDEXES FOR FASTER SEARCHES
        cmd = """ CREATE INDEX idx_group ON events (groupName)"""
        cursor.execute(cmd,)
        
        cmd = """ CREATE INDEX idx_uid ON events (uid) """
        cursor.execute(cmd,)
            
    return None


def create_first_song_played():

    

    tableName = 'events'
    with sqlite3.connect(sqlitePath) as conn:
        cursor = conn.cursor()

        cmd = """
        CREATE TABLE first_song AS
        SELECT uid, min(DATE) AS date
        FROM events
        WHERE event = 'song_played'
        GROUP BY uid
        """

        cursor.execute(cmd,)


def create_subscribed_table():


    tableName = 'events'
    with sqlite3.connect(sqlitePath) as conn:
        cursor = conn.cursor()

        cmd = """
        CREATE TABLE  subscribed  AS
        SELECT uid,DATE
        FROM events
        WHERE event = 'subscribed'
        """

        cursor.execute(cmd,)

    return None


def create_unsubscription_table():


    tableName = 'events'
    with sqlite3.connect(sqlitePath) as conn:
        cursor = conn.cursor()

        cmd = """
        CREATE TABLE  unsubscribed  AS
        SELECT uid,DATE
        FROM events
        WHERE event = 'unsubscribed'
        """

        cursor.execute(cmd,)

    return None

    
def create_account_creation_table():


    
    with sqlite3.connect(sqlitePath) as conn:
        cursor = conn.cursor()

        cmd = """
        CREATE TABLE  account  AS
        SELECT uid,DATE
        FROM events
        WHERE event = 'account created'
        """

        cursor.execute(cmd,)

    return None
