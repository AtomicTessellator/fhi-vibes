import sqlite3

def convertToColType(inputVal, lenIdVals=1):
    if(type(inputVal) is tuple):
        return [ [ inputVal ] ]
    elif(type(inputVal) is list):
        if(len(inputVal) > 0):
            if(type(inputVal[0]) is list):
                return inputVal
            elif(type(inputVal[0]) is tuple):
                newInputs = []
                for ii in inputVal:
                    newInputs.append([ii])
                return newInputs
            else:
                raise TypeError("For adding/updating a row column a tuple of name/value must be given")
        else:
            for ii in range(lenIdVals):
                inputVal.append([])
            return inputVal
    else:
        raise TypeError("For adding/updating a row column a tuple of name/value must be given")

def convert2KeyVal(inputVal):
    if(type(inputVal) == list):
        return inputVal
    else:
        return [inputVal]

def adapt_list(ll):
    frmt = ""
    for l in ll:
        frmt += f"{l};{};"
    return frmt[:-1].encode()

def convert_list(bs):
    bsSplit = bs.split(b";")
    ll = []


class MateirailsDatabase:
    def __init__(self, filename, adapt_reg=[], con_reg=[]):
        for adap in adapt_reg:
            sqlite3.register_adapter(adap[0], adap[1] )
        for con in con_reg:
            sqlite3.register_converter(con[0], con[1] )
        self.fileNameDB_ = filename
        self.connDB_ = sqlite3.connect(self.fileNameDB_, detect_types=sqlite3.PARSE_DECLTYPES )
        self.cursorDB_ = self.connDB_.cursor()
        self.tableNameList_ = self.getTableList()

    def addAdapter2Register(self, obj, adapt_fxn):
        sqlite3.register_adapter(obj, adapt_fxn )

    def addConverter2Register(self, sqlLiteName, convert_fxn):
        sqlite3.register_converter(sqlLiteName, convert_fxn )

    def getTableList(self):
        self.cursorDB_.execute("SELECT name FROM sqlite_master WHERE type='table';")
        tables = self.cursorDB_.fetchall()
        return tables

    def createTable(self, tableName, cols):
        makeTableString = "CREATE TABLE IF NOT EXISTS " + tableName + " ("
        for col in cols:
            makeTableString += col[0] + " " + col[1] + ", "
        makeTableString = makeTableString[:-2] + ")"
        self.cursorDB_.execute(makeTableString)
        self.tableNameList_.append(tableName)
        return

    def selectRows(self, tableName, idCol, idVals):
        idVals = convert2KeyVal(idVals)
        rows = []
        for idVal in idVals:
            self.cursorDB_.execute('SELECT * FROM ' + tableName + ' WHERE ' + idCol + '=?', (idVal,) )
            rows.append(self.cursorDB_.fetchall())
        return rows

    def selectCols(self, tableName, colIDs):
        colIDs = convert2KeyVal(colIDs)
        selectString = "SELECT " + colIDs[0]
        for ii in range(1, len(colIDs)):
            selectString += ", " + colIDs[ii]
        selectString += " FROM " + tableName
        self.cursorDB_.execute(selectString)
        cols = self.cursorDB_.fetchall()
        for tt in range(len(cols)):
            if(len(cols[tt]) == 1):
                cols[tt] = cols[tt][0]
        return cols

    def selectCell(self, tableName, idCol, idVal, colName):
        self.cursorDB_.execute("SELECT " + colName + " FROM " + tableName + " WHERE " + idCol + "=?", (idVal,) )
        return self.cursorDB_.fetchone()[0]

    def insertRows(self, tableName, idCol, idVals, colNamesVals):
        idVals = convert2KeyVal(idVals)
        colNamesVals = convertToColType(colNamesVals, len(idVals))
        if(len(idVals) != len(colNamesVals)):
            raise VlaueError("length of idVals and colNamesVals, must be the same")
        for ii in range(len(colNamesVals)):
            try:
                colNameStr = " (" + idCol
                colValTuple = (idVals[ii],)
                for col in colNamesVals[ii]:
                    colNameStr += ", " + col[0]
                    colValTuple += (col[1],)
                colNameStr += ") VALUES (?"
                for jj in range(len(colNamesVals[ii])):
                    colNameStr += ", ?"
                colNameStr += ")"
                insertString = "INSERT INTO " + tableName + colNameStr
                self.cursorDB_.execute(insertString, colValTuple)
            except sqlite3.IntegrityError:
                print('ERROR: ID already exists in PRIMARY KEY column ' + idCol + ": " + idVals[ii])
        return

    def insertCols(self, tableName, colNamesTypes):
        for col in colNamesTypes:
            self.cursorDB_.execute("ALTER TABLE " + tableName + " ADD COLUMN '" + col[0] + "' " + col[1])
        return

    def updateRow(self, tableName, idCol, idVals, colNamesVals):
        idVals = convert2KeyVal(idVals)
        colNamesVals = convertToColType(colNamesVals, len(idVals))
        if(len(idVals) != len(colNamesVals)):
            raise VlaueError("length of idVals and colNamesVals, must be the same")
        for ii in range(len(colNamesVals)):
            colNameStr = ""
            colValTuple = ()
            for jj in range(len(colNamesVals[ii])):
                colNameStr += colNamesVals[ii][jj][0] + "=? , "
                colValTuple += (colNamesVals[ii][jj][1],)
            colNameStr = colNameStr[:-3]
            colValTuple += (idVals[ii],)
            self.cursorDB_.execute("UPDATE " + tableName + " SET " + colNameStr + " WHERE " + idCol + "=?", colValTuple )
        return

    def deleteRows(self, tableName, idCol, idVals):
        idVals = convert2KeyVal(idVals)
        for val in idVals:
            self.cursorDB_.execute("DELETE FROM " + tableName + " WHERE " + idCol + "=" + val)
        return

    def dropTable(self, tableName):
        self.cursorDB_.execute("DROP TABLE IF EXISTS " + tableName)

    def close(self):
        self.cursorDB_.close()
        self.connDB_.close()

    def execute(self, execStr):
        self.cursorDB_.execute(execStr)