import sqlite3

# Had an Idea to maintain a database instead of processing the data every time

conn = sqlite3.connect('HWG.db')

cursor = conn.cursor()

conn.commit()
conn.close()