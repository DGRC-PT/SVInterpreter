#!/usr/bin/python3
encoding="UTF-8"

def index2():
        print('content-type:text/html; charset=utf-8 \n\n')
        print ('''
<html>
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<link href="https://fonts.googleapis.com/css?family=Arial" rel="stylesheet">
<head>
<meta charset="UTF-8">
<title>SVInterpreter</title>
<link rel="stylesheet" href="/outputs/test.css">
<script src="/outputs/JS_script.js"></script>
<link href="/outputs/nouislider.min.css" rel="stylesheet">
</head>
<body>
<div class="container">
<br><br/>
<center><h1>SVInterpreter - Search Results</h1></center>
<hgh><center><p><a3>The retrieved data from each breakpoint is compiled in a table that can be vizualized on the browser or downloaded in xlsx.</a3></p></center></hgh>
<center><hgh><button onclick="goBack()">New search</button></center></hgh>
<center><b>The results will appear shortly....</b></center>
<script>
function goBack() {
  window.location.replace("../SVInterpreter.py");
}
</script>
</div>
</body>
</html>
''')

index2()
