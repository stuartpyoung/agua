dojo.provide("plugins.core.Common");

/* 	CLASS SUMMARY: PROVIDE COMMONLY USED METHODS FOR ALL CLASSES.

	ALSO PROVIDE LOW-LEVEL METHODS THAT ACCOMPLISH GENERIC TASKS WHICH

	ARE WRAPPED AROUND BY CONTEXT-SPECIFIC METHODS
*/

// COMMON METHODS FOR ALL PLUGINS
dojo.declare( "plugins.core.Common",
    null,
{
	//////}

// HASH OF LOADED CSS FILES
loadedCssFiles : null,

constructor : function () {
},


// *************
// HASHARRAY METHODS
// *************
getObjectByKeyValue : function (hasharray, key, value) {
// RETRIEVE AN OBJECT FROM AN ARRAY, IDENTIFED BY A KEY:VALUE PAIR



	if ( hasharray == null )	return;

	for ( var i = 0; i < hasharray.length; i++ )
	{
		if ( hasharray[i][key] == value )	return hasharray[i];
	}

	return;
},

_addObjectToArray : function (array, object, requiredKeys) {
// ADD AN OBJECT TO AN ARRAY, CHECK REQUIRED KEY VALUES ARE NOT NULL


	if ( array == null )
	{
		array = new Array;
	}

	var notDefined = this.notDefined(object, requiredKeys);
	if ( notDefined.length > 0 )
	{
		return false;
	}

	array.push(object);

	return true;
},

_removeObjectFromArray : function (array, object, keys) {
// REMOVE AN OBJECT FROM AN ARRAY, IDENTIFY OBJECT USING SPECIFIED KEY VALUES


	// DEBUG splitace.pl
	var arrayCopy = dojo.clone(array);

	var index = this._getIndexInArray(array, object, keys);

	if ( index == null )	return false;

	array.splice(index, 1);

	return true;
},

_removeArrayFromArray : function(hasharray, removeThis, uniqueKeys) {
// REMOVE AN ARRAY OF OBJECTS FROM A LARGER ARRAY CONTAINING IT,
// IDENTIFYING OBJECTS USING THE SPECIFIED KEY VALUES


	var success = true;
	for ( var j = 0; j < removeThis.length; j++ )
	{

		var notDefined = this.notDefined (removeThis[j], uniqueKeys);
		if ( notDefined.length > 0 )
		{
			success = false;
			continue;
		}

		var removeSuccess = this._removeObjectFromArray(hasharray, removeThis[j], uniqueKeys);

		////////// ***** DEBUG ONLY *****
		////////if ( removeSuccess == false )
		////////{				
		////////	success = false;
		////////}
	}

	return success;
},

_getIndexInArray : function (hasharray, object, keys) {
// GET THE INDEX OF AN OBJECT IN AN ARRAY, IDENTIFY OBJECT USING SPECIFIED KEY VALUES


	if ( hasharray == null )
	{
		return null;
	}

	for ( var i = 0; i < hasharray.length; i++ )
	{
		var arrayObject = hasharray[i];
		var identified = true;
		for ( var j = 0; j < keys.length; j++ )
		{
			if ( arrayObject[keys[j]] != object[keys[j]] )
			{
				identified = false;
				break;
			}
		}

		if ( identified == true )
		{
			return i;
		}
	}

	return null;
},

_objectInArray : function (array, object, keys) {
// RETURN true IF AN OBJECT ALREADY BELONGS TO A GROUP


	if ( array == null )	return false;
	if ( object == null )	return false;

	var index = this._getIndexInArray(array, object, keys);

	if ( index == null )	return false;

	return true;
},


// *************
// two-D ARRAY METHODS
// *************
_addArrayToArray : function (twoDArray, array, requiredKeys) {
// ADD AN ARRAY TO A TWO-D ARRAY, CHECK REQUIRED KEY VALUES ARE NOT NULL


	if ( twoDArray == null )
	{
		twoDArray = new Array;
	}

	var notDefined = this.notDefined(array, requiredKeys);
	if ( notDefined.length > 0 )
	{
		return false;
	}

	//var satisfied = true;
	//for ( var i in requiredKeys )
	//{
	//	if ( array[requiredKeys[i]] == null )
	//	{
	//		satisfied = false;
	//		break;
	//	}
	//}
	//
	//if ( satisfied == false )
	//{
	//	return false;
	//}
	//

	twoDArray.push(array);

	return true;
},


// ******************	
// DATA METHODS
// ******************	
notDefined : function (hasharray, keys) {
// RETURN AN ARRAY OF THE UNDEFINED FIELDS IN A HASH

	if ( hasharray == null )	return;
	var notDefined = new Array;
	for ( var i = 0; i < keys.length; i++ )
	{
		if ( hasharray[keys[i]] == null )	notDefined.push(keys[i]);
	}

	return notDefined;
},

sortHasharrayByOrder : function (hasharray, order) {
// SORT A HASHARRAY BY THE GIVEN ORDER OF KEYS, EXCLUDING ENTRIES
// THAT DO NOT HAVE VALUES FOR ANY OF THE GIVEN KEYS


	var orderedHasharray = new Array;

	for ( var i = 0; i < order.length; i++ )
	{
		var orderedType = order[i];

		for ( var j = 0; j < hasharray.length; j++ )
		{
			var applicationHash = hasharray[j];
			var applicationType;
			for ( var applicationType in applicationHash )
			{
				if ( applicationType == orderedType )
				{
					orderedHasharray.push(applicationHash);
					break;
				}
			}
		}
	}

	return orderedHasharray;		
},

sortNoCase : function (array) {
// DO A NON-CASE SPECIFIC SORT OF AN ARRAY


	return array.sort( function (a,b)
		{
			return a.toUpperCase() == b.toUpperCase() ?
			(a < b ? -1 : a > b) : (a.toUpperCase() < b.toUpperCase() ? -1 : a.toUpperCase() > b.toUpperCase());
		}
	);
},

hasharrayKeyToArray : function (hasharray, key) {
// RETURN AN ARRAY CONTAINING THE SPECIFIED key VALUE
// IN EACH HASH IN A HASHARRAY

	var outputArray = new Array;
	dojo.forEach(hasharray, function(entry) {
		if ( entry[key] != null )
		{
			outputArray.push(entry);
		}
	});

	return outputArray;
},

hashkeysToArray : function (hash) {
// RETURN AN ARRAY CONTAINING ONLY THE KEYS OF A HASH

	var array = new Array;
	for ( var key in hash )
	{
		if ( key != null && key != '' )   array.push(key);
	}

	return array;
},

filterHasharray : function (hasharray, key) {
// RETURN A HASHARRAY CONTAINING THE SPECIFIED key VALUE
// IN EACH HASH IN A HASHARRAY

	var outputHasharray = new Array;
	dojo.forEach(hasharray, function(entry) {
		if ( entry[key] != null )
		{
			outputHasharray.push(entry);
		}
	});

	return outputHasharray;
},

filterByKeyValues : function (hasharray, keyarray, valuearray ) {
// RETURN AN ARRAY OF OBJECTS THAT ALL POSSESS GIVEN KEY VALUE
// NB: THIS SPLICES OUT THE ENTRIES FROM THE ACTUAL INPUT ARRAY
// REFERENCE I.E., THE PASSED ARRAY WILL SHRINK IN SIZE

	if ( hasharray == null )	return;
	if ( keyarray == null )	return;
	if ( valuearray == null )	return;

	for ( var i = 0; i < hasharray.length; i++ )
	{
		var isMatched = true;
		for ( var j = 0; j < keyarray.length; j++ )
		{
			if ( hasharray[i][keyarray[j]] != valuearray[j] )
			{
				isMatched = false;
				break;
			}
		}
		if ( isMatched == false )	{ hasharray.splice(i, 1);	i--; }
	}

	return hasharray;
},

sortHasharrayByKeys : function (hashArray, keys) {

	if ( hashArray == null )	return;
	if ( keys == null )	return;
	if ( ! typeof hashArray == "ARRAY" || hashArray == null ) return;
	if ( ! typeof keys == "ARRAY" || keys == null ) return;

	return hashArray.sort(function (a,b) {



        var result = 1;
        for ( var i = 0; i < keys.length; i++ )
        {
            var key = keys[i];
            if ( a[key] == null || b[key] == null )
            {
                continue;
            }

            if ( a[key].toUpperCase && b[key].toUpperCase )
            {

                var aString = a[key].toUpperCase(); 
                var bString = b[key].toUpperCase(); 

                result =  aString == bString ?
                ( a[key] < b[key] ? -1 : (a[key] > b[key] ? 1 : 0) )
                    : ( aString < bString ? -1 : (aString > bString ? 1 : 0) );

                if ( result != 0 )    break;
            }
            else
            {
                result = a[key] < b[key] ? -1 :
                         ( a[key] > b[key] ? 1 : 0 );

                if ( result != 0 )    break;
            }  
        }

        return result;
    });
},


sortHasharray : function (hashArray, key) {
// SORT AN ARRAY OF HASHES BY A SPECIFIED HASH FIELD
// NB: IF THE FIELD IS NULL OR EMPTY IN AN ARRAY ENTRY
// IT WILL BE DISCARDED.

	if ( hashArray == null )	return;
	if ( key == null )	return;
	if ( ! typeof hashArray == "ARRAY" ) return;

	//dojo.toJson(hashArray);


//	if ( hashArray == null )
//	{
//		return;
//	}
//	
//	if ( key == null )
//	{
//		return;
//	}


	return hashArray.sort(function (a,b) {

		if ( a[key] == null )
		{
			return;
		}
		if ( b[key] == null )
		{
			return;
		}

		return a[key].toUpperCase() == b[key].toUpperCase() ?
			(a[key] < b[key] ? -1
			: a[key] > b[key])
			: (a[key].toUpperCase() < b[key].toUpperCase() ? -1
			: a[key].toUpperCase() > b[key].toUpperCase());
		}
	);
},

sortNumericHasharray : function (hashArray, key) {
// SORT AN ARRAY OF HASHES BY A SPECIFIED *NUMERIC* HASH FIELD
// NB: IF THE FIELD IS NULL OR EMPTY IN AN ARRAY ENTRY
// IT WILL BE DISCARDED.

//		if ( hashArray == null )	return;
//		if ( key == null )	return;

	if ( hashArray == null )
	{
		return;
	}

	if ( key == null )
	{
		return;
	}

	// REMOVE NON-NUMERIC ENTRIES FOR SORT KEY
	for ( var i = 0; i < hashArray.length; i++ )
	{
		if ( parseInt(hashArray[i]) == "NaN" )
		{
			hashArray.splice(i, 1);
			i--;
		}
	}

	return hashArray.sort(function (a,b) {
		return parseInt(a[key]) == parseInt(b[key]) ?
			(  parseInt(a[key]) < parseInt(b[key]) ? -1
			: parseInt(a[key]) > parseInt(b[key])  )
			: (  parseInt(a[key]) < parseInt(b[key]) ? -1
			: parseInt(a[key]) > parseInt(b[key]) );
		}
	);

},

sortTwoDArray : function (twoDArray, index) {
// SORT AN ARRAY OF HASHES BY A SPECIFIED HASH FIELD
// NB: IF THE FIELD IS NULL OR EMPTY IN AN ARRAY ENTRY
// IT WILL BE DISCARDED.


//	if ( twoDArray == null )
//	{
//		return;
//	}
//	
//	if ( index == null )
//	{
//		return;
//	}
	if ( twoDArray == null )	return;
	if ( index == null )	return;

	return twoDArray.sort(function (a,b) {

		if ( a[index] == null )
		{
			return;
		}
		if ( b[index] == null )
		{
			return;
		}

		return a[index].toUpperCase() == b[index].toUpperCase() ?
			(a[index] < b[index] ? -1
			: a[index] > b[index])
			: (a[index].toUpperCase() < b[index].toUpperCase() ? -1
			: a[index].toUpperCase() > b[index].toUpperCase());
		}
	);

},

hasharrayKeyToArray : function (hasharray, key) {
// PUT A SPECIFIED HASHARRAY KEY INTO AN ARRAY

if ( hasharray == null )	return;
	if ( key == null )	return;

	var array = new Array;
	for ( var i = 0; i < hasharray.length; i++ )
	{
		if ( hasharray[i][key] != null 
			&& hasharray[i][key] != '' )
					array.push(hasharray[i][key]);
	}

	return array;
},

getEntry : function (hasharray, keyarray, valuearray ) {
// RETURN AN ENTRY IN A HASH ARRAY IDENTIFIED BY ITS UNIQUE KEYS

	if ( hasharray == null )	return;
	if ( keyarray == null )	return;
	if ( valuearray == null )	return;

	for ( var i = 0; i < hasharray.length; i++ )
	{
		var isMatched = true;
		for ( var j = 0; j < keyarray.length; j++ )
		{
			if ( hasharray[i][keyarray[j]] != valuearray[j] )
			{
				isMatched = false;
				break;
			}
		}

		if ( isMatched == true )
		{
			return hasharray[i];
		}
	}

	return null;
},


uniqueValues : function(array) {

	if ( array.length == 1 )	return array;

	array = array.sort();		
	for ( var i = 1; i < array.length; i++ )
	{
		if ( array[i-1] == array[i] )
		{
			array.splice(i, 1);
			i--;

		}
	}

	return array;
},



// *************
// COMBOBOX METHODS
// *************
setUsernameCombo : function () {
//	POPULATE COMBOBOX AND SET SELECTED ITEM
//	INPUTS: Agua.sharedprojects DATA OBJECT
//	OUTPUTS:	ARRAY OF USERNAMES IN COMBO BOX, ONCLICK CALL TO setSharedProjectCombo

	//console.dir(Agua);

	var usernames = Agua.getSharedUsernames();

	// RETURN IF projects NOT DEFINED
	if ( usernames == null || usernames.length == 0 )
	{
		return;
	}

	// DO data FOR store
	var data = {identifier: "name", items: []};
	for ( var i in usernames )
	{
		data.items[i] = { name: usernames[i]	};
	}

	// CREATE store
	var store = new dojo.data.ItemFileReadStore( {	data: data	} );

	// ADD STORE TO USERNAMES COMBO
	this.usernameCombo.store = store;	

	// START UP AND SET VALUE
	this.usernameCombo.startup();
	this.usernameCombo.set('value', usernames[0]);			
},

setSharedProjectCombo : function (username, projectName, workflowName) {
//	POPULATE COMBOBOX AND SET SELECTED ITEM
//	INPUTS: USERNAME, OPTIONAL PROJECT NAME AND WORKFLOW NAME
//	OUTPUTS: ARRAY OF USERNAMES IN COMBO BOX, ONCLICK CALL TO setSharedWorkflowCombo


	var projects = Agua.getSharedProjectsByUsername(username);
	if ( projects == null )
	{
		return;
	}

	var projectNames = this.hasharrayKeyToArray(projects, "project");
	projectNames = this.uniqueValues(projectNames);

	// RETURN IF projects NOT DEFINED
	if ( projectNames == null || projectNames.length == 0 )
	{
		return;
	}

	// DO data FOR store
	var data = {identifier: "name", items: []};
	for ( var i in projectNames )
	{
		data.items[i] = { name: projectNames[i]	};
	}

	// CREATE store
	var store = new dojo.data.ItemFileReadStore( {	data: data	} );

	// ADD STORE TO USERNAMES COMBO
	this.projectCombo.store = store;	

	// START UP AND SET VALUE
	this.projectCombo.startup();
	this.projectCombo.set('value', projectNames[0]);	
},

setSharedWorkflowCombo : function (username, projectName, workflowName) {
//	POPULATE COMBOBOX AND SET SELECTED ITEM
//	INPUTS: USERNAME, OPTIONAL PROJECT NAME AND WORKFLOW NAME
//	OUTPUTS: ARRAY OF USERNAMES IN COMBO BOX, ONCLICK CALL TO setSharedWorkflowCombo


	if ( projectName == null )	projectName = this.projectCombo.get('value');

	var workflows = Agua.getSharedWorkflowsByProject(username, projectName);
	if ( workflows == null )
	{
		return;
	}

	var workflowNames = this.hasharrayKeyToArray(workflows, "workflow");
	workflowNames = this.uniqueValues(workflowNames);

	// RETURN IF workflows NOT DEFINED
	if ( workflowNames == null || workflowNames.length == 0 )
	{
		return;
	}

	// DO data FOR store
	var data = {identifier: "name", items: []};
	for ( var i in workflowNames )
	{
		data.items[i] = { name: workflowNames[i]	};
	}

	// CREATE store
	var store = new dojo.data.ItemFileReadStore( {	data: data	} );

	// ADD STORE TO USERNAMES COMBO
	this.workflowCombo.store = store;	

	// START UP AND SET VALUE
	this.workflowCombo.startup();
	this.workflowCombo.set('value', workflowNames[0]);
},

setProjectCombo : function (project, workflow) {
//	INPUT: (OPTIONAL) project, workflow NAMES
//	OUTPUT:	POPULATE COMBOBOX AND SET SELECTED ITEM


	var projectNames = Agua.getProjectNames();

	// RETURN IF projects NOT DEFINED
	if ( ! projectNames )
	{
		return;
	}

	// SET PROJECT IF NOT DEFINED TO FIRST ENTRY IN projects
	if ( project == null || ! project)	project = projectNames[0];

	// DO DATA ARRAY
	var data = {identifier: "name", items: []};
	for ( var i in projectNames )
	{
		data.items[i] = { name: projectNames[i]	};
	}

	// CREATE store
	var store = new dojo.data.ItemFileReadStore( {	data: data	} );

	//// GET PROJECT COMBO WIDGET
	var projectCombo = this.projectCombo;
	if ( projectCombo == null )
	{
		return;
	}

	projectCombo.store = store;	

	// START UP AND SET VALUE
	//projectCombo.startup();
	projectCombo.setValue(project);			
},

setWorkflowCombo : function (project, workflow) {
// SET THE workflow COMBOBOX


	if ( project == null || ! project )
	{
		return;
	}

	// CREATE THE DATA FOR A STORE		
	var workflows = Agua.getWorkflowNamesByProject(project);


	// RETURN IF workflows NOT DEFINED
	if ( ! workflows )
	{
		return;
	}		

	// DO data FOR store
	var data = {identifier: "name", items: []};
	for ( var i in workflows )
	{
		data.items[i] = { name: workflows[i]	};
	}

	// CREATE store
	var store = new dojo.data.ItemFileReadStore( { data: data } );

	// GET WORKFLOW COMBO
	var workflowCombo = this.workflowCombo;
	if ( workflowCombo == null )
	{
		return;
	}

	workflowCombo.store = store;

	// START UP COMBO AND SET SELECTED VALUE TO FIRST ENTRY IN workflows IF NOT DEFINED 
	if ( workflow == null || ! workflow )	workflow = workflows[0];

	workflowCombo.startup();
	workflowCombo.set('value', workflow);			
},

setReportCombo : function (project, workflow, report) {
// SET THE report COMBOBOX


	if ( project == null || ! project )
	{
		return;
	}
	if ( workflow == null || ! workflow )
	{
		return;
	}

	var reports = Agua.getReportsByProjectWorkflow(project, workflow);
	if ( reports == null )	reports = [];

	var reportNames = this.hasharrayKeyToArray(reports, "name");

	// DO data FOR store
	var data = {identifier: "name", items: []};
	for ( var i in reports )
	{
		data.items[i] = { name: reportNames[i]	};
	}

	// CREATE store
	// http://docs.dojocampus.org/dojo/data/ItemFileWriteStore
	var store = new dojo.data.ItemFileReadStore( { data: data } );

	// GET WORKFLOW COMBO
	var reportCombo = this.reportCombo;
	if ( reportCombo == null )
	{
		return;
	}

	reportCombo.store = store;

	// GET USER INPUT WORKFLOW
	var snpReport = this;

	// START UP COMBO (?? NEEDED ??)
	reportCombo.startup();
	reportCombo.set('value', report);			
},

setViewCombo : function (projectName, viewName) {
// SET THE view COMBOBOX

	// SANITY CHECK
	if ( this.viewCombo == null )	return;

	if ( projectName == null || ! projectName )
	{
		return;
	}

	var views = Agua.viewNames(projectName);
	views = views.sort();

	// RETURN IF views NOT DEFINED
	if ( ! views || views.length == 0 )
	{
		return;
		//Agua.addView({ project: projectName, name: "View1" });
		//views = Agua.viewNames(projectName);
	}		

	// SET view IF NOT DEFINED TO FIRST ENTRY IN views
	if ( viewName == null || ! viewName)
	{
		viewName = views[0];
	}

	// DO data FOR store
	var data = {identifier: "name", items: []};
	for ( var i in views )
	{
		data.items[i] = { name: views[i]	};
	}

	// CREATE store
	// http://docs.dojocampus.org/dojo/data/ItemFileWriteStore
	var store = new dojo.data.ItemFileReadStore( { data: data } );

	this.viewCombo.store = store;

	// START UP COMBO (?? NEEDED ??)
	this.viewCombo.startup();
	this.viewCombo.set('value', viewName);			
},

// *************
// HOUSEKEEPING METHODS
// *************
loadCSS : function (cssFiles) {
// LOAD EITHER this.cssFiles OR A SUPPLIED ccsFiles ARRAY OF FILES ARGUMENT


	if ( cssFiles == null )
	{
		cssFiles = this.cssFiles;
	}

	// LOAD CSS
	for ( var i in cssFiles )
	{

		this.loadCSSFile(cssFiles[i]);
	}
},

loadCSSFile : function (cssFile) {
// LOAD A CSS FILE IF NOT ALREADY LOADED, REGISTER IN this.loadedCssFiles


	if ( Agua.loadedCssFiles == null || ! Agua.loadedCssFiles )
	{
		Agua.loadedCssFiles = new Object;
	}

	if ( ! Agua.loadedCssFiles[cssFile] )
	{

		var cssNode = document.createElement('link');
		cssNode.type = 'text/css';
		cssNode.rel = 'stylesheet';
		cssNode.href = cssFile;
		document.getElementsByTagName("head")[0].appendChild(cssNode);

		Agua.loadedCssFiles[cssFile] = 1;
	}
	else
	{
	}


	return Agua.loadedCssFiles;
},


systemVariables : function (value) {
// PARSE THE SYSTEM VARIABLES: %TEXT% VALUES


	value = value.replace(/%project%/g, this.project);
	value = value.replace(/%workflow%/g, this.workflow);

	return value;
},

insertTextBreaks : function (text) {
// INSERT INVISIBLE UNICODE CHARACTER &#8203; AT INTERVALS OF
// LENGTH width IN THE TEXT


	// SET INSERT CHARACTERS
	var insert = "\n";

	// FIRST, REMOVE ANY EXISTING INVISIBLE LINE BREAKS IN THE TEXT
	text = this.removeTextBreaks(text);

	// SECOND, INSERT A "&#8203;" CHARACTER AT REGULAR INTERVALS
	var insertedText = '';
	var offset = 0;
	while ( offset < text.length )
	{
		var temp = text.substring(offset, offset + this.textBreakWidth);
		offset += this.textBreakWidth;

		insertedText += temp;
		//insertedText += "&#8203;"
		insertedText += insert;
	}

	return insertedText;
},

removeTextBreaks : function (text) {
// REMOVE ANY LINE BREAK CHARACTERS IN THE TEXT


	text = text.replace(/\n/g, '');


	text = text.replace(/\n/g, '');

	return text;
},

jsonSafe : function (string, mode) {
// CONVERT FROM JSON SAFE TO ORDINARY TEXT OR THE REVERSE


	// SANITY CHECKS
	if ( string == null || ! string )
	{
		return '';
	}
	if ( string == true || string == false )
	{
		return string;
	}

	// CLEAN UP WHITESPACE
	string = String(string).replace(/\s+$/g,'');
	string = String(string).replace(/^\s+/g,'');
	string = String(string).replace(/"/g,"'");

	var specialChars = [
		[ '&quot;',	"'" ],	
		//[ '&quot;',	'"' ],	
		[ '&#35;', '#' ],	
		[ '&#36;', '$' ],	
		[ '&#37;', '%' ],	
		[ '&amp;', '&' ],	
		//[ '&#39;', "'" ],	
		[ '&#40;', '(' ],	
		[ '&#41;', ')' ],	
		[ '&frasl;', '/' ],	
		[ '&#91;', '\[' ],	
		[ '&#92;', '\\' ],	
		[ '&#93;', '\]' ],	
		[ '&#96;', '`' ],	
		[ '&#123;', '\{' ],	
		[ '&#124;', '|' ],	
		[ '&#125;', '\}' ]	
	];

	for ( var i = 0; i < specialChars.length; i++)
	{
		if ( mode == 'toJson' )
		{
			var sRegExInput = new RegExp(specialChars[i][1], "g");
			return string.replace(sRegExInput, specialChars[i][0]);
		}
		else if ( mode == 'fromJson' )
		{
			var sRegExInput = new RegExp(specialChars[i][0], "g");
			return string.replace(sRegExInput, specialChars[i][1]);
		}
	}

	return string;
},

autoHeight : function(textarea) {

	var rows = parseInt(textarea.rows);

	while ( textarea.rows * textarea.cols < textarea.value.length )
	{
		textarea.rows = ++rows;
	}

	// REMOVE ONE LINE TO FIT SNUGLY
	//textarea.rows--; 

////    // THIS ALSO WORKS OKAY
////    
////    var height = event.target.scrollHeight;
////    var rows = parseInt((height / 15) - 2);
////
////
//////<textarea id="value" class="autosize" rows="3" cols="10">    
////    
////    
////    event.target.setAttribute('rows', rows);
////    
////
},

firstLetterUpperCase : function (word) {
	if ( word == null || word == '' ) return word;
	if ( word.substring == null )	return null;
	return word.substring(0,1).toUpperCase() + word.substring(1);
},


cleanEnds : function (words) {
	if ( words == null )	return '';

	words = words.replace(/^[\s\n]+/, '');
	words = words.replace(/[\s\n]+$/, '');

	return words;
},


cleanWord : function (word) {

	if ( word == null )	return '';

	return word.replace(/[\s\n]+/, '');
},

clearValue : function (textarea, value) {
	if ( textarea == null )	return;
	if ( textarea.value == null )	return;
	if ( textarea.value.match(value) )
	{
		textarea.value = '';
	}
},





quit : function (message) {
	return;
},


randomiseUrl : function (url) {
	url += "?dojo.preventCache=1266358799763";
	url += Math.floor(Math.random()*1000000000000);

	return url;
},


doPut : function (inputs) {
	var callback = function (){};
	if ( inputs.callback != null )	callback = inputs.callback;
	var url = inputs.url;
	url += "?";
	url += Math.floor(Math.random()*100000);
	var query = inputs.query;
	var timeout = inputs.timeout ? inputs.timeout : 30000;
	//var handleAs = inputs.handleAs ? inputs.handleAs : "json-comment-filtered";
	var handleAs = inputs.handleAs ? inputs.handleAs : "text";
	var sync = inputs.sync ? inputs.sync : false;

	// SEND TO SERVER
	dojo.xhrPut(
		{
			url: url,
			contentType: "text",
			preventCache : true,
			sync: sync,
			handleAs: handleAs,
			putData: dojo.toJson(query),
			timeout: timeout,
			load: function(response, ioArgs) {
				callback();
			},
			error: function(response, ioArgs) {
				return response;
			}
		}
	);	
},


killPopup : function (combo) {
	var popupId = "widget_" + combo.id + "_dropdown";
	var popup = dojo.byId(popupId);
	if ( popup != null )
	{
		var popupWidget = dijit.byNode(popup.childNodes[0]);
		dijit.popup.close(popupWidget);
	}	
}




});

