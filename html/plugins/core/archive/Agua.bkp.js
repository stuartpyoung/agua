dojo.provide( "plugins.core.Agua");

/*

	PURPOSE:

		1. PROVIDE INTERFACE WITH Agua DATA OBJECT REPRESENTATION

			OF THE DATA MODEL ON THE REMOTE SERVER

		2. PROVIDE METHODS TO CHANGE/INTERROGATE THE DATA OBJECT

		3. CALLS TO REMOTE SERVER TO REFLECT CHANGES ARE MOSTLY THE

			RESPONSIBILITY OF THE OBJECT USING THE Agua CLASS

*/

// INHERITS
dojo.require("plugins.core.Common");
dojo.require("plugins.core.Updater");
dojo.require("dijit.Toolbar");
dojo.require("dijit.layout.TabContainer");

// TOASTER
dojo.require("dojox.widget.Toaster");
dojo.require("dijit.Tooltip");

dojo.require("dojox.widget.Standby");

dojo.declare( "plugins.core.Agua",
	[ dijit._Widget, dijit._Templated, plugins.core.Common ],
{
//Path to the template of this widget. 
templatePath: dojo.moduleUrl("plugins", "core/templates/agua.html"),	

// CSS files
cssFiles : [ "dojo-1.5.0/dojox/layout/resources/ExpandoPane.css",
			 "plugins/core/css/agua.css" ],

// Calls dijit._Templated.widgetsInTemplate
widgetsInTemplate : true,


////}


// PLUGINS TO LOAD
// NB: ORDER IS IMPORTANT FOR CORRECT LAYOUT
pluginsList : [
	//"plugins.login.Controller"
	//,
	//"plugins.admin.Controller"
	//,
	//"plugins.project.Controller"
	//,
	"plugins.workflow.Controller"
	//,
	//"plugins.report.Controller"
	//,
	//"plugins.view.Controller"
	//,
	//"plugins.help.Controller"
],

////{

// CONTROLLERS
controllers : new Object(),

// DIV FOR PRELOAD SCREEN
splashNode : null,

// DIV TO DISPLAY PRELOAD MESSAGE BEFORE MODULES ARE LOADED
messageNode : null,

// PLUGIN MANAGER LOADS THE PLUGINS
pluginManager: null,

// COOKIES CONTAINS STORED USER ID AND SESSION ID
cookies : new Object,

// CONTAINS ALL LOADED CSS FILES
css : new Object,

// WEB URLs
cgiUrl : null,
htmlUrl : null,

////}

// CONSTRUCTOR
constructor : function(args) {


	// SET DIVS FOR PRELOAD SCREEN AND PRELOAD MESSAGE
	this.splashNode = args.splashNode;
	this.messageNode = args.messageNode;
	this.cgiUrl = args.cgiUrl;
	this.htmlUrl = args.htmlUrl;

},

postCreate: function() {

	this.startup();
},

addToolbarButton: function (label) {
// ADD MODULE BUTTON TO TOOLBAR


	if ( this.toolbar == null )
	{
		return;
	}

	var button = new dijit.form.Button({
		label: label,
		showLabel: true,
		iconClass: "dijitEditorIcon dijitEditorIcon" + label
	});
	this.toolbar.addChild(button);

	return button;
},


startup : function () {
// CHECK IF DEBUGGING AND LOAD PLUGINS


	this.loadCSS();




	// ATTACH THIS TEMPLATE TO attachPoint DIV ON HTML PAGE
	var attachPoint = dojo.byId("attachPoint");
	attachPoint.appendChild(this.containerNode);

	// LOAD FIRST BUTTON IN TOOLBAR
	this.aguaButton = this.addToolbarButton("Agua");

	// SET BUTTON LISTENER
	var listener = dojo.connect(this.aguaButton, "onClick", this, "reload");

	// SET UP THE ELEMENT OBJECTS AND THEIR VALUE FUNCTIONS
	this.inherited(arguments);


	// INITIALISE ELECTIVE UPDATER
	this.updater = new plugins.core.Updater();




	//<div id ="standby10" target="overlayTarget10" dojoType="dojox.widget.Standby"></div>


	//// SET STANDBY
	//this.standby = new dojox.widget.Standby(
	//	{
	//		target: this.controls
	//		//,
	//		//onClick: "reload",
	//		//id : "standby"
	//	}
	//);
	//
	//this.standby.show();
	////console.dir(this.standby);
	//

},



reload : function () {
// RELOAD AGUA


	var url = window.location;
	window.open(location, '_blank', 'toolbar=1,location=0,directories=0,status=0,menubar=1,scrollbars=1,resizable=1,navigation=0'); 

	//window.location.reload();
},



loadPlugins : function () {
	dojo.require("plugins.core.PluginManager");

	// CHECK IF DEBUGGING
	var isDebug = this.debug();



	// LOAD PLUGINS
	this.showPreloadMessage("Loading plugins");
	this.pluginManager = new plugins.core.PluginManager({
		parentWidget : this,
		pluginsList : this.pluginsList
	})

	var url = window.location.href;
	if ( isDebug == true && ! url.match(/(.+?)\?([^\?]+),([^,]+)$/) )
	{
		var baseUrl = url.match(/(.+?)\?([^\?]+),([^,]+)$/)[1];
		window.location.href = baseUrl + "?test,999";
	}

//return;
//
	//if ( this.loginController == null )
	//{
	//	this.workflowController = new plugins.workflow.Controller();
	//	
	//	// DEBUG
	//	this.workflowController.createTab();
	//}
},


debug : function () {
// IF DEBUGGING, SET USERNAME AND SESSIONID FROM URL AND LOAD DATA


	// GET USERNAME FROM URL IF PRESENT
	var url = window.location.href;
	if ( ! url.match(/(.+?)\?([^\?]+),([^,]+)$/) )
	{
		return false;
	}

	var username = url.match(/(.+?)\?([^\?]+),([^,]+)$/)[2];
	var sessionId = url.match(/(.+?)\?([^\?]+),([^,]+)$/)[3];

	// SET AGUA COOKIE
	this.cookie('username', username);
	this.cookie('sessionId', sessionId);


	// GET DATA FROM SERVER
	this.showPreloadMessage("Loading data");
	this.getData();

	// REMOVE plugins.login.Controller FROM this.pluginsList
	for ( var i = 0; i < this.pluginsList.length; i++ )
	{
		if ( this.pluginsList[i] == "plugins.login.Controller" )
		{
			this.pluginsList.splice(i, 1);
		}
	}


	return true;
},

toastMessage : function (message, type, duration) {

	if ( duration == null )	duration = 1000;
	if ( type != null && ( type != "warning" && type != "error" ) )
	{
		return;
	}

	var topic = "toastTopic";
	dojo.publish(topic, [ {
		message: message,
		type: type,
		duration: duration
	}]);
},


startPlugins : function () {

	// GET DATA FROM SERVER
	this.showPreloadMessage("Loading data");
	this.getData();

	this.adminController = new plugins.admin.Controller();
	this.adminController.createTab();

	this.projectController = new plugins.project.Controller();
	this.projectController.createTab();

	this.workflowController = new plugins.workflow.Controller();
	this.workflowController.createTab();

	//this.reportController = new plugins.report.Controller();
	//this.reportController.createTab();
	//
	//Agua.viewController = new plugins.view.Controller();
	//Agua.viewController.createTab();

	this.helpController = new plugins.help.Controller();
	this.helpController.createTab();
},

// DATA METHODS
getData : function() {
// LOAD ALL AGUA DATA FOR THIS USER, INCLUDING SHARED PROJECT DATA


/*
this.headings = {
	leftPane : [ "Access" ],
	middlePane : [],
	rightPane : []
};

 this.access = [
	{
	  'owner' : 'syoung',
	  'worldcopy' : '1',
	  'groupwrite' : '1',
	  'worldwrite' : '1',
	  'groupname' : 'mihg',
	  'groupcopy' : '1',
	  'groupview' : '1',
	  'worldview' : '1'
	},
	{
	  'owner' : 'syoung',
	  'worldcopy' : '1',
	  'groupwrite' : '1',
	  'worldwrite' : '1',
	  'groupname' : 'snp',
	  'groupcopy' : '1',
	  'groupview' : '1',
	  'worldview' : '1'
	}
];

return;
*/

	// GET URL 
	var url = this.cgiUrl + "workflow.cgi";

	// GENERATE QUERY JSON FOR THIS WORKFLOW IN THIS PROJECT
	var query = new Object;
	query.username = this.cookie('username');
	query.sessionId = this.cookie('sessionId');
	query.mode = "getData";

	var aguaObject = this;
	dojo.xhrPut(
		{

			url: url,
			putData: dojo.toJson(query),
			//handleAs: "json",
			handleAs: "json-comment-optional",
			sync: true,
			load: function(response, ioArgs) {
				if ( response.error )
				{

				}
				else
				{
					for ( var key in response )
					{
						// ARRAY OF PROJECT HASHES
						aguaObject[key] = response[key];
					}
				}
			},
			error: function(response, ioArgs) {
			}
		}
	);
},

// ADMIN METHODS
getHeadings : function () {

	return dojo.clone(this.headings)
},

getAccess : function () {

	return dojo.clone(this.access)
},

// FILE METHODS
getFileExists : function (stageParameterObject, booleanValue) {
// GET THE BOOLEAN fileExists VALUE FOR A STAGE PARAMETER


	if ( booleanValue != null )
	{
		return null;
	}

	return this._fileExists(stageParameterObject, booleanValue);
},


setFileExists : function (stageParameterObject, booleanValue) {
// SET THE BOOLEAN fileExists VALUE FOR A STAGE PARAMETER


	if ( booleanValue == null )
	{
		return null;
	}

	return this._fileExists(stageParameterObject, booleanValue);
},

_fileExists : function (stageParameterObject, booleanValue) {
// RETURN THE fileExists BOOLEAN FOR A STAGE PARAMETER
// OR SET IT IF A VALUE IS SUPPLIED


	//var filtered = dojo.clone(this.stageparameters);
	//var keys = ["appname"];
	//var values = ["image2eland.pl"];
	//filtered = this.filterByKeyValues(filtered, keys, values);


	var uniqueKeys = ['username', 'project', 'workflow', 'appname', 'appnumber', 'name'];
	var valueArray = new Array;
	for ( var i = 0; i < uniqueKeys.length; i++ )
	{
		valueArray.push(stageParameterObject[uniqueKeys[i]]);
	}
	var stageParameter = this.getEntry(this.stageparameters, uniqueKeys, valueArray);
	if ( stageParameter == null )
	{
		return null;
	}

	if ( booleanValue != null )
	{

		// SET exists BOOLEAN VALUE
		stageParameter.exists = booleanValue;

//console.dir(stageParameter);

		var success = this._removeStageParameter(stageParameter);
		if ( success == false )
		{
			return null;
		}



		success = this._addStageParameter(stageParameter);			
		if ( success == false )
		{
			return null;
		}

		return true;
	}
	else
	{
		return stageParameter.exists;
	}
},

// SOURCE METHODS
getSources : function () {
// RETURN A SORTED COPY OF this.sources


	var sources = dojo.clone(this.sources)
	return this.sortHasharray(sources, "name");
},


isSource : function (sourceObject) {
// RETURN TRUE IF SOURCE NAME ALREADY EXISTS


	var sources = this.getSources();
	if ( sources == null )	return false;

	return this._objectInArray(sources, sourceObject, ['name']);
},

addSource : function (sourceObject) {
// ADD A SOURCE OBJECT TO THE this.sources ARRAY


	return this._addObjectToArray(this.sources, sourceObject, [ 'name', 'description', 'location' ]);
},


removeSource : function (sourceObject) {
// REMOVE A SOURCE OBJECT FROM this.sources AND this.groupmembers


	var success = this._removeObjectFromArray(this.sources, sourceObject, ['name']);
	if ( success == null || ! success )
	{
	}

	// REMOVING SOURCE FROM groups TABLE

	// ADD USERNAME AND TYPE TO SOURCE OBJECT
	sourceObject.username = Agua.cookie('username');
	sourceObject.type = "source";

	var requiredKeys = [ "username", "name", "type"];
	return this._removeObjectFromArray(this.groupmembers, sourceObject, requiredKeys);
},

isGroupSource : function (groupName, sourceObject) {
// RETURN true IF A SOURCE ALREADY BELONGS TO A GROUP


	var groupSources = this.getGroupSources();
	if ( groupSources == null )	return false;

	groupSources = this.filterByKeyValues(groupSources, ["groupname"], [groupName]);

	return this._objectInArray(groupSources, sourceObject, ['name']);
},

addSourceToGroup : function (groupName, sourceObject) {
// ADD A SOURCE OBJECT TO A GROUP ARRAY IF IT DOESN'T EXIST THERE ALREADY 


	if ( this.isGroupSource(groupName, sourceObject) == true )
	{
		return false;
	}

	var groups = this.getGroups();
	if ( groups == null )	return null;
	var group = this.getObjectByKeyValue(groups, "name", groupName);

	sourceObject.username = group.username;
	sourceObject.groupname = groupName;
	sourceObject.groupdesc = group.description;
	sourceObject.type = "source";

	var requiredKeys = [ "username", "groupname", "name", "type"];
	return this._addObjectToArray(this.groupmembers, sourceObject, requiredKeys);
},

removeSourceFromGroup : function (groupName, sourceObject) {
// REMOVE A SOURCE OBJECT FROM A GROUP ARRAY, IDENTIFY OBJECT BY 'name' KEY VALUE


	var groups = this.getGroups();
	if ( groups == null )	return null;
	var group = this.getObjectByKeyValue(groups, "name", groupName);

	sourceObject.username = group.username;
	sourceObject.groupname = groupName;
	sourceObject.groupdesc = group.description;
	sourceObject.type = "source";

	var requiredKeys = [ "username", "groupname", "name", "type"];
	return this._removeObjectFromArray(this.groupmembers, sourceObject, requiredKeys);
},



getSourcesByGroup : function (groupName) {
// RETURN THE ARRAY OF SOURCES THAT BELONG TO A GROUP

	var groupSources = this.getGroupSources();
	if ( groupSources == null )	return null;

	var keyArray = ["groupname"];
	var valueArray = [groupName];
	return this.filterByKeyValues(groupSources, keyArray, valueArray);
},


// AWS METHODS
getAws : function (username) {
// RETURN ENTRY FOR username IN this.aws

	return dojo.clone(this.aws);
},


setAws : function (aws) {
// RETURN ENTRY FOR username IN this.aws

	if ( aws == null )
	{
		return;
	}
	if ( aws.amazonuserid == null )
	{
		return;
	}

	this.aws = aws;

	return aws;
},


// USER METHODS
getUser : function (username) {
// RETURN ENTRY FOR username IN this.users

	var users = this.getUsers();

	var index = this._getIndexInArray(users, [username], [0]);
	if ( index != null )
	{
		return this.users[index];
	}

	return null;
},

getUsers : function () {
// RETURN A SORTED COPY OF this.users

	var users = dojo.clone(this.users);

	return this.sortHasharray(users, 0);
},

addUser : function (userObject) {
// ADD A USER OBJECT TO THE this.users ARRAY


	// ARRAY FORMAT:
	// userArray[0]: ["aabate","a","abate","aabate@med.miami.edu",""]
	var userArray = new Array;
	userArray[0] = userObject.username;
	userArray[1] = userObject.firstname || '';
	userArray[2] = userObject.lastname || '';
	userArray[3] = userObject.email || '';

	// NEW USER MUST HAVE username AND email
	var removed = this._removeArrayFromArray(this.users, userArray, [0]);

	var added = this._addArrayToArray(this.users, userArray, [ 0, 3 ]);
	if ( added == false )
	{
		return false;
	}

	this.users = this.sortHasharray(this.users, 0);

	return true;
},


isUser : function (userObject) {
// ADD A USER OBJECT TO THE this.users ARRAY


	// ARRAY FORMAT:
	// userArray[0]: ["aabate","a","abate","aabate@med.miami.edu",""]
	var userArray = new Array;
	userArray[0] = userObject.username;
	userArray[1] = userObject.firstname || '';
	userArray[2] = userObject.lastname || '';
	userArray[3] = userObject.email || '';

	if ( this._getIndexInArray(this.users, userArray, [0]) )	{	return 1;	}

	return 0;
},

removeUser : function (userObject) {
// REMOVE A USER OBJECT FROM THE this.users ARRAY


	// ARRAY FORMAT:
	// userArray[0]: ["aabate","a","abate","aabate@med.miami.edu",""]
	var userArray = new Array;
	userArray[0] = userObject.username;
	userArray[1] = userObject.firstname || '';
	userArray[2] = userObject.lastname || '';
	userArray[3] = userObject.email || '';

	// DELETED USER MUST HAVE username
	return this._removeArrayFromArray(this.users, userArray, [0]);
},


isGroupUser : function (groupName, userObject) {
// RETURN true IF A USER ALREADY BELONGS TO A GROUP


	var groupUsers = this.getGroupUsers();
	if ( groupUsers == null )	return false;

	groupUsers = this.filterByKeyValues(groupUsers, ["groupname"], [groupName]);

	userObject.name = userObject.username;

	return this._objectInArray(groupUsers, userObject, ['name']);
},


addUserToGroup : function (groupName, userObject) {
// ADD A USER OBJECT TO A GROUP ARRAY IF IT DOESN'T EXIST THERE ALREADY 

	if ( this.isGroupUser(groupName, userObject) == true )
	{
		return false;
	}

	var groups = this.getGroups();
	if ( groups == null )	return null;
	var group = this.getObjectByKeyValue(groups, "name", groupName);

	userObject.username = group.username;
	userObject.groupname = groupName;
	userObject.groupdesc = group.description;
	userObject.type = "user";

	var requiredKeys = [ "username", "groupname", "name", "type"];
	return this._addObjectToArray(this.groupmembers, userObject, requiredKeys);
},


removeUserFromGroup : function (groupName, userObject) {
// REMOVE A USER OBJECT FROM A GROUP ARRAY, IDENTIFY OBJECT BY 'name' KEY VALUE


	var groups = this.getGroups();
	if ( groups == null )	return null;
	var group = this.getObjectByKeyValue(groups, "name", groupName);

	userObject.owner = group.username;
	userObject.groupname = groupName;
	userObject.groupdesc = group.description;
	userObject.type = "user";

	var requiredKeys = [ "username", "groupname", "name", "type"];
	return this._removeObjectFromArray(this.groupmembers, userObject, requiredKeys);
},




// PROJECT METHODS
getProjects : function () {
// RETURN A SORTED COPY OF this.projects


	var projects = dojo.clone(this.projects)
	return this.sortHasharray(projects, "name");
},

getProjectNames : function (projects) {
// RETURN AN ARRAY OF ALL PROJECT NAMES IN this.projects


	if ( projects == null )	projects = this.projects;
	if ( projects == null )	return;

	return this.hasharrayKeyToArray(projects, "name");
},


addProject : function (projectObject) {
// ADD A PROJECT OBJECT TO THE this.projects ARRAY


	if ( this.addingProject == true )	return;
	this.addingProject == true;

	//var projectObject = new Object;
	//projectObject.name = projectName;
	// ADD USERNAME TO PROJECT OBJECT
	projectObject.username = Agua.cookie('username');

	var added = this._addObjectToArray(this.projects, projectObject, [ 'name' ]);
	if ( added == false )
	{
		this.addingProject == false;
		return;
	}

	this.addingProject == false;

	//// COMMIT CHANGES IN REMOTE DATABASE
	//var url = Agua.cgiUrl + "workflow.cgi";
	//
	//// SET QUERY
	//var query = new Object;
	//query.data = projectObject;
	//query.username = Agua.cookie('username');
	//query.sessionId = Agua.cookie('sessionId');
	//query.mode = "addProject";
	//
	//this.doPut({ url: url, query: query, sync: false });
},

copyProject : function (sourceUser, sourceProject, targetUser, targetProject, copyFiles) {
// ADD AN EMPTY NEW WORKFLOW OBJECT TO A PROJECT OBJECT


	if ( this.isProject(targetProject) == true )
	{
		return;
	}

	// ADD WORKFLOW TO TOP OF WORKFLOW ARRAY IN PROJECT
	var projectObject = new Object;
	projectObject.name = targetProject;
	var keys = ["name"];
	var copied = this._addObjectToArray(this.projects, projectObject, keys);
	if ( copied == false )
	{
		return;
	}

	// COMMIT CHANGES TO REMOTE DATABASE
	var url = Agua.cgiUrl + "agua";
	var query = new Object;
	query.sourceuser = sourceUser;
	query.targetuser = targetUser;
	query.sourceproject = sourceProject;
	query.targetproject = targetProject;
	query.copyfiles = copyFiles;
	query.username = Agua.cookie('username');
	query.sessionId = Agua.cookie('sessionId');
	query.mode = "copyProject";

	this.doPut({ url: url, query: query, sync: false });
},


removeProject : function (projectObject) {
// REMOVE A PROJECT OBJECT FROM: this.projects, this.workflows, this.groupmembers
// this.stages AND this.stageparameters


	// RETURN IF name IS NULL
	if ( projectObject.name == null )
	{
		return;
	}

	// SET ADDITIONAL FIELDS
	projectObject.project = projectObject.name;
	projectObject.owner = Agua.cookie('username');
	projectObject.type = "project";

	// REMOVE PROJECT FROM this.projects
	var success = this._removeObjectFromArray(this.projects, projectObject, ['name']);

	// REMOVE FROM this.workflows
	var success = this._removeObjectFromArray(this.workflows, projectObject, ['project']);

	// REMOVE FROM this.groupmembers 
	var keys = ["owner", "name", "type"];
	this._removeObjectFromArray(this.groupmembers, projectObject, keys);

	// REMOVE FROM this.stages AND this.stageparameters
	var keys = ["project"];
	this._removeObjectFromArray(this.stages, projectObject, keys);

	// REMOVE FROM this.stageparameters
	var keys = ["project"];
	this._removeObjectFromArray(this.stageparameters, projectObject, keys);

	//// COMMIT CHANGES IN REMOTE DATABASE
	//var url = Agua.cgiUrl + "workflow.cgi";
	//var query = projectObject;
	//query.username = Agua.cookie('username');
	//query.sessionId = Agua.cookie('sessionId');
	//query.mode = "removeProject";
	//
	//this.doPut({ url: url, query: query, sync: false });
},


isProject : function (projectName) {
// RETURN true IF A PROJECT EXISTS IN this.projects


	if ( this.projects == null )	return false;
	for ( var i in this.projects )
	{
		var project = this.projects[i];

		if ( project.name.toLowerCase() == projectName.toLowerCase() )
		{
			return true;
		}
	}

	return false;
},


isGroupProject : function (groupName, projectObject) {
// RETURN true IF A PROJECT BELONGS TO THE SPECIFIED GROUP


	var groupProjects = this.getGroupProjects();
	if ( groupProjects == null )	return false;

	groupProjects = this.filterByKeyValues(groupProjects, ["groupname"], [groupName]);

	return this._objectInArray(groupProjects, projectObject, ['groupname', 'name']);
},


addProjectToGroup : function (groupName, projectObject) {
// ADD A PROJECT OBJECT TO A GROUP ARRAY IF IT DOESN'T EXIST THERE ALREADY 


	if ( this.isGroupProject(groupName, projectObject) == true )
	{
		return false;
	}


	projectObject.owner = projectObject.username;
	projectObject.groupname = groupName;
	projectObject.type = "project";


	// ********* DEBUG **********
	// REMOVE THIS LATER		
	// ********* DEBUG **********

	var groups = this.getGroups();
	if ( groups == null )	return null;
	var group = this.getObjectByKeyValue(groups, "name", groupName);
	projectObject.groupdesc = group.description;


	var requiredKeys = ["owner", "groupname", "name", "type"];
	return this._addObjectToArray(this.groupmembers, projectObject, requiredKeys);
},


removeProjectFromGroup : function (groupName, projectObject) {
// REMOVE A PROJECT OBJECT FROM A GROUP ARRAY, IDENTIFY OBJECT BY 'name' KEY VALUE


	projectObject.owner = projectObject.username;
	projectObject.groupname = groupName;
	projectObject.type = "project";



	// REMOVE THIS LATER		
	var groups = this.getGroups();
	if ( groups == null )	return null;
	var group = this.getObjectByKeyValue(groups, "name", groupName);
	projectObject.groupdesc = group.description;



	var requiredKeys = [ "owner", "groupname", "name", "type"];
	return this._removeObjectFromArray(this.groupmembers, projectObject, requiredKeys);
},


getProjectsByGroup : function (groupName) {
// RETURN THE ARRAY OF PROJECTS THAT BELONG TO A GROUP

	var groupProjects = this.getGroupProjects();
	if ( groupProjects == null )	return null;

	var keyArray = ["groupname"];
	var valueArray = [groupName];
	var projects = this.filterByKeyValues(groupProjects, keyArray, valueArray);

	return this.sortHasharray(projects, "name");
},

getGroupsByProject : function (projectName) {
// RETURN THE ARRAY OF PROJECTS THAT BELONG TO A GROUP

	var groupProjects = this.getGroupProjects();
	if ( groupProjects == null )	return null;

	var keyArray = ["type", "name"];
	var valueArray = ["project", projectName];
	return this.filterByKeyValues(groupProjects, keyArray, valueArray);
},


// 	WORKFLOW METHODS
getWorkflows : function () {
// RETURN A SORTED COPY OF this.workflows


	var workflows = dojo.clone(this.workflows)
	return this.sortHasharray(workflows, "name");
},

getWorkflowNamesByProject : function (projectName) {
// RETURN AN ARRAY NAMES OF WORKFLOWS IN this.workflows WITH THE SPECIFIED PROJECT


	if ( this.workflows == null )	return;

	var workflows = dojo.clone(this.workflows);
	var keyArray = ["project"];
	var valueArray = [projectName];
	workflows = this.filterByKeyValues(workflows, keyArray, valueArray);

	workflows = this.sortHasharray(workflows, "name");

	return this.hasharrayKeyToArray(workflows, "name");
},

getWorkflow : function (projectName, workflowName) {
// RETURN TRUE IF WORKFLOW NAME IS FOUND IN PROJECT IN this.workflows


	var result = this._getWorkflow(projectName, workflowName);
	if ( result == null || ! result.length || result.length == 0 )
	{
		return null;
	}

	var workflowObject = result[0];	

	return workflowObject;
},


_getWorkflow : function (projectName, workflowName) {
// RETURN WORKFLOW IN this.workflows IDENTIFIED BY PROJECT AND WORKFLOW NAMES


	if ( this.isProject(projectName) == false )
	{
		return false;
	}

	// GET ALL WORKFLOWS
	var workflows = this.getWorkflows();
	if ( workflows == null ) 
	{
		return false;
	}

	// CHECK FOR OUR PROJECT AND WORKFLOW NAME AMONG WORKFLOWS
	var keyArray = ["project", "name"];
	var valueArray = [projectName, workflowName];
	return this.filterByKeyValues(workflows, keyArray, valueArray);
},


isWorkflow : function (projectName, workflowName) {
// RETURN TRUE IF WORKFLOW NAME IS FOUND IN PROJECT IN this.workflows


	var result = this._getWorkflow(projectName, workflowName);
	if ( result == null || ! result.length || result.length == 0 )
	{
		return false;
	}

	return true;		
},


addWorkflow : function (projectName, workflowName) {
// ADD AN EMPTY NEW WORKFLOW OBJECT TO A PROJECT OBJECT


	if ( this.isWorkflow(projectName, workflowName)== true )
	{
		return;
	}

	// ADD WORKFLOW TO TOP OF WORKFLOW ARRAY IN PROJECT
	var workflowObject = new Object;
	workflowObject.name = workflowName;
	workflowObject.project = projectName;

	var keys = ["project", "name"];
	var added = this._addObjectToArray(this.workflows, workflowObject, keys);
	if ( added == false )
	{
		return;
	}

	// COMMIT CHANGES IN REMOTE DATABASE
	var url = Agua.cgiUrl + "workflow.cgi";

	// SET QUERY
	var query = new Object;
	query.name = workflowName;
	query.project = projectName;
	query.username = Agua.cookie('username');
	query.sessionId = Agua.cookie('sessionId');
	query.mode = "addWorkflow";

	this.doPut({ url: url, query: query, sync: false });
},

copyWorkflow : function (sourceUser, sourceProject, sourceWorkflow, targetUser, targetProject, targetWorkflow, copyFiles) {
// ADD AN EMPTY NEW WORKFLOW OBJECT TO A PROJECT OBJECT


	if ( this.isWorkflow(targetProject, targetWorkflow)== true )
	{
		return;
	}

	this._copyWorkflow(sourceUser, sourceProject, sourceWorkflow, targetUser, targetProject, targetWorkflow, copyFiles);
	// ADD STAGES, STAGEPARAMETERS, REPORTS AND VIEWS



	// COMMIT CHANGES TO REMOTE DATABASE
	var url = Agua.cgiUrl + "workflow.cgi";
	var query = new Object;
	query.sourceuser = sourceUser;
	query.targetuser = targetUser;
	query.sourceworkflow = sourceWorkflow;
	query.sourceproject = sourceProject;
	query.targetworkflow = targetWorkflow;
	query.targetproject = targetProject;
	query.copyfiles = copyFiles;
	query.username = Agua.cookie('username');
	query.sessionId = Agua.cookie('sessionId');
	query.mode = "copyWorkflow";

	this.doPut({ url: url, query: query, sync: false });
},

_copyWorkflow : function (sourceUser, sourceProject, sourceWorkflow, targetUser, targetProject, targetWorkflow, copyFiles) {

	// COPY WORKFLOW 
	var workflowObject = new Object;
	workflowObject.name = targetWorkflow;
	workflowObject.project = targetProject;
	var keys = ["project", "name"];
	var copied = this._addObjectToArray(this.workflows, workflowObject, keys);
	if ( copied == false )
	{
		return;
	}

	// COPY STAGES AND STAGE PARAMETERS
	var stages;
	if ( sourceUser != targetUser )
		stages = this.getSharedStagesByWorkflow(sourceUser, sourceProject, sourceWorkflow);
	else stages = this.getStagesByWorkflow(sourceProject, sourceWorkflow);
	for ( var i = 0; i < stages.length; i++ )
	{
		this._addStage(stages[i]);

		// ADD STAGE PARAMETERS
		var stageparams;
		if ( sourceUser != targetUser )
			stageparams = this.getSharedStageParameters(stages[i]);
		else stageparams = this.getStageParameters(stages[i]);
		for ( var i = 0; i < stageparams.length; i++ )
		{
			this._addStageParameter(stageparams[i]);
		}
	}

	// COPY VIEWS
	var views = this.getSharedViews({ username: sourceUser, project: sourceProject } );
	for ( var i = 0; i < views.length; i++ )
	{
		this._addStage(views[i]);
	}


},


removeWorkflow : function (workflowObject) {
// REMOVE A WORKFLOW FROM this.workflows, this.stages AND this.stageparameters


	if ( workflowObject.name == null )	return;
	if ( workflowObject.project == null )	return;

	// REMOVE FROM this.workflows
	var uniqueKeys = ['project', 'name'];
	var result = this._removeObjectFromArray(this.workflows, workflowObject, uniqueKeys);


	// REMOVE FROM this.stages AND this.stageparameters
	workflowObject.workflow = workflowObject.name;
	var keys = ["project", "workflow"];
	this._removeObjectFromArray(this.stages, workflowObject, keys);
	this._removeObjectFromArray(this.stageparameters, workflowObject, keys);


	// COMMIT CHANGES IN REMOTE DATABASE
	var url = Agua.cgiUrl + "workflow.cgi";
	var query = new Object;
	query.project = workflowObject.project;
	query.name = workflowObject.name;
	query.username = Agua.cookie('username');
	query.sessionId = Agua.cookie('sessionId');
	query.mode = "removeWorkflow";

	this.doPut({ url: url, query: query, sync: false });
},


renameWorkflow : function (workflowObject, newName) {
// RENAME A WORKFLOW FROM this.workflows, this.stages AND this.stageparameters


	if ( workflowObject.name == null )	return;
	if ( workflowObject.project == null )	return;
	if ( newName == null )	return;

	// SAVE OLD NAME
	var oldName = workflowObject.name;

	// RENAME FROM this.workflows
	var uniqueKeys = ['project', 'name'];
	var result = this._removeObjectFromArray(this.workflows, workflowObject, uniqueKeys);

	// RENAME FROM this.stages AND this.stageparameters
	workflowObject.workflow = workflowObject.name;

	// REMOVE WORKFLOW OBJECT AND ADD IT BACK WITH NEW NAME
	var keys = ["project", "workflow"];
	this._removeObjectFromArray(this.stages, workflowObject, keys);
	workflowObject.name = newName;
	this._addObjectToArray(this.stageparameters, workflowObject, keys);

	// COMMIT CHANGES IN REMOTE DATABASE
	var url = Agua.cgiUrl + "workflow.cgi";
	var query = new Object;
	query.project = workflowObject.project;
	query.name = oldName;
	query.newname = newName;
	query.username = Agua.cookie('username');
	query.sessionId = Agua.cookie('sessionId');
	query.mode = "renameWorkflow";

	this.doPut({ url: url, query: query, sync: false });
},



// GROUP METHODS
getGroups : function () {
// RETURN THE this.groups ARRAY FOR THIS USER


	return dojo.clone(this.groups);
},

addGroup : function (groupObject) {
// ADD A GROUP OBJECT TO THE this.groups ARRAY


	// DO THE ADD
	var result = this._addObjectToArray(this.groups, groupObject, [ 'name' ]);

	// ADD GROUP TO groupMembers HASH
	if ( result == true ) this.groupmembers[groupObject.name] = new Object;

	// SORT GROUPS ARRAY
	if ( result == true ) this.sortGroups();

	// RETURN TRUE OR FALSE
	return result;
},


removeGroup : function (groupObject) {
// REMOVE A GROUP OBJECT FROM THE this.groups ARRAY AND THE this.groupmembers ARRAY


	var result = this._removeObjectFromArray(this.groups, groupObject, ['name']);

	// REMOVE GROUP FROM groupMembers HASH
	if ( result == true ) delete this.groupmembers[groupObject.name];

	// SORT GROUPS ARRAY
	if ( result == true ) this.sortGroups();

	return result;
},


isGroup : function (groupName) {
// RETURN true IF A GROUP EXISTS IN this.groups


	if ( this.groups == null )	return false;
	for ( var i in this.groups )
	{
		var group = this.groups[i];
		if ( group.name.toLowerCase() == groupName.toLowerCase() )
		{
			return true;
		}
	}

	return false;
},




sortGroups : function () {
// DO NON-CASE SPECIFIC SORT OF this.groups BY GROUP NAME

	if ( this.groups == null )
	{
		return;
	}

	// DO CASE-INSENSITIVE SORT AND PUT SORTED ARRAY IN this.groups
	this.groups = this.sortHasharray(this.groups, "name");

	if ( this.groups != null )	return true;
	return false;
},


getGroupNames : function () {
// PARSE NAMES OF ALL GROUPS IN this.groups INTO AN ARRAY

	if ( this.groups == null ) return [];

	var groupNames = new Array;
	for ( var i in this.groups  )
	{
		groupNames.push(this.groups[i].name);
	}

	return groupNames;
},


getGroupMembers : function (memberType) {
// PARSE this.groups ENTRIES INTO HASH OF ARRAYS groupName: [ source1, source2 ]


	if ( this.groupmembers == null )
	{
		return;
	}

	var groupMembers = dojo.clone(this.groupmembers);
	var keyArray = ["type"];
	var valueArray = [memberType];
	return this.filterByKeyValues(groupMembers, keyArray, valueArray);
},


getGroupMembersHash : function (memberType) {
// PARSE this.groups ENTRIES INTO HASH OF ARRAYS { groupName: [ source1, source2 ] }


	if ( this.groups == null )
	{
		return;
	}

	var groupMembers = dojo.clone(this.groupmembers);
	for ( var groupName in groupMembers )
	{
		for ( var j = 0; j < groupMembers[groupName].length; j++ )
		{
			if ( groupMembers[groupName][j].type != memberType )
			{
				groupMembers[groupName].splice(j,1);
				j--;
			}
		} 
	}

	return groupMembers;
},

getGroupSources : function () {
// GET ALL SOURCE MEMBERS OF this.groupmembers

	return this.getGroupMembers('source');
},


getGroupUsers : function () {
// GET ALL USER MEMBERS OF this.groupmembers

	return this.getGroupMembers('user');
},


getGroupProjects : function () {
// GET ALL PROJECT MEMBERS OF this.groupmembers

	return this.getGroupMembers('project');
},








// SHARED METHODS
getSharedSources : function () {
// RETURN A COPY OF this.sharedsources
	if ( this.sharedsources != null )
		return dojo.clone(this.sharedsources);

	else return [];
},

getSharedUsernames : function() {
// RETURN AN ARRAY OF ALL OF THE NAMES OF USERS WHO HAVE SHARED PROJECTS
// WITH THE LOGGED ON USER


	// PUT THE USER NAMES IN THIS ARRAY
	if ( this.sharedprojects == null || this.sharedprojects.length == 0 )
	{
	}

	var usernames = new Array;
	var alreadySeen = new Object;
	for ( var i = 0; i < this.sharedprojects.length; i++ )
	{
		if ( ! alreadySeen[this.sharedprojects[i].username] )
		{
			alreadySeen[this.sharedprojects[i].username] = 1;
			usernames.push(this.sharedprojects[i].username);
		}
	}

	return usernames;
},


getSharedProjectsByUsername : function(username) {
/// RETURN AN ARRAY OF ALL OF THE NAMES OF USERS WHO HAVE SHARED PROJECTS
// WITH THE LOGGED ON USER


	var sharedProjects = dojo.clone(this.sharedprojects);

	if ( sharedProjects == null )	return [];

	sharedProjects = this.filterByKeyValues(sharedProjects, ["username"], [username]);


	return sharedProjects;
},


getSharedStagesByUsername : function(username) {
// RETURN AN ARRAY OF STAGES SHARED WITH THIS USER 

	if ( this.sharedstages == null )	return [];

	var sharedStages = dojo.clone(this.sharedstages[username]);
	if ( sharedStages == null )	return [];

	return sharedStages;
},


getSharedWorkflowsByProject : function(username, project) {
// RETURN AN ARRAY OF ALL OF THE NAMES OF USERS WHO HAVE SHARED PROJECTS
// WITH THE LOGGED ON USER


	var sharedStages = this.getSharedStagesByUsername(username);
	sharedStages = this.filterByKeyValues(sharedStages, ["project"], [project]);


	return sharedStages;
},



getSharedStagesByWorkflow : function(username, project, workflow) {
// RETURN AN ARRAY OF ALL OF THE NAMES OF USERS WHO HAVE SHARED PROJECTS
// WITH THE LOGGED ON USER


	var sharedStages = this.getSharedWorkflowsByProject(username, project);

	sharedStages = this.filterByKeyValues(sharedStages, ["workflow"], [workflow]);


	return sharedStages;
},

getSharedStageParameters : function (stageObject) {
getSharedParametersByAppname : function (appname, owner) {
// RETURN AN ARRAY OF PARAMETERS FOR THE GIVEN APPLICATION


	if ( appname == null )	return null;	
	var parameters = new Array;
	var params = dojo.clone(this.parameters);
	dojo.forEach(params, function(parameter)
	{
		if ( parameter.appname == appname ) parameters.push(parameter);
	});

	return parameters;
},

// RETURN AN ARRAY OF STAGE PARAMETER HASHARRAYS FOR A STAGE


	var keys = ["username", "project", "workflow", "name", "number"];
	var notDefined = this.notDefined (stageObject, keys);
	if ( notDefined.length != 0 )
	{
		return;
	}

	var stageParameters = dojo.clone(this.sharedstageparameters[stageObject.username]);

	// ADD appname AND appnumber FOR STAGE PARAMETER IDENTIFICATION
	stageObject.appname = stageObject.name;
	stageObject.appnumber = stageObject.number;

	var keyArray = ["username", "project", "workflow", "appname", "appnumber"];
	var valueArray = [stageObject.username, stageObject.project, stageObject.workflow, stageObject.name, stageObject.number];
	stageParameters = this.filterByKeyValues(stageParameters, keyArray, valueArray);


	return stageParameters;
},


getSharedViews : function (viewObject) {
// RETURN AN ARRAY OF VIEWS IN THE SHARED PROJECT

	var keys = ["project", "username"];
	var notDefined = this.notDefined (viewObject, keys);
	if ( notDefined.length != 0 )
	{
		return;
	}

	var views = dojo.clone(this.sharedviews[viewObject.username]);

	var keyArray = ["project", "username"];
	var valueArray = [viewObject.username, viewObject.project, viewObject.workflow, viewObject.name, viewObject.number];
	views = this.filterByKeyValues(views, keyArray, valueArray);


	return views;
},

// APP METHODS
getAppTypes : function (appGroup) {
// GET SORTED LIST OF ALL APP TYPES

	var apps = this.getApps(appGroup);

	var typesHash = new Object;
	for ( var i = 0; i < apps.length; i++ )
	{
		typesHash[apps[i].type] = 1;
	}

	var types = this.hashkeysToArray(typesHash)
	types = this.sortNoCase(types);

	return types;
},

sortApps : function (appGroup) {	
// PARSE NAMES OF ALL GROUPS IN this.apps[appGroup] INTO AN ARRAY

	if ( this.apps == null ){
		return;
	}

	if ( this.apps[appGroup] == null )
	{
		return;
	}

	// DO CASE-INSENSITIVE SORT AND PUT SORTED ARRAY IN this.apps[appGroup]
	this.apps[appGroup] = this.sortHasharray(this.apps[appGroup], "name");

	if ( this.apps[appGroup] != null )	return true;
	return false;
},

hasApps : function (appGroup) {
// ADD AN APP OBJECT TO THE this.apps[appGroup] OR this.apps.common ARRAY


	var hasApps = this.apps[appGroup] != null ? true : false;
	return hasApps;
},

getApps : function (appGroup) {
// RETURN THE this.apps[appGroup] ARRAY FOR THIS USER


	if ( this.apps == null ) 	return [];
	if ( this.apps[appGroup] == null ) 	return [];

	return dojo.clone(this.apps[appGroup]);
},

getAppsByType : function (type, appGroup) {
// ADD AN APP OBJECT TO THE this.apps[appGroup] ARRAY


	var apps = this.getApps(appGroup);
	var keyArray = ["type"];
	var valueArray = [type];
	apps = this.filterByKeyValues(apps, keyArray, valueArray);

	return apps;
},

addApp : function (appObject, appGroup) {
// ADD AN APP OBJECT TO THE this.apps[appGroup] ARRAY


	// DO THE ADD
	var result = this._addObjectToArray(this.apps[appGroup], appObject, [ 'name' ]);

	// SORT this.appNames ARRAY
	if ( result == true ) this.sortApps();

	// RETURN TRUE OR FALSE
	return result;
},

removeApp : function (appObject, appGroup) {
// REMOVE AN APP OBJECT FROM THE this.apps[appGroup] ARRAY 


	var result = this._removeObjectFromArray(this.apps[appGroup], appObject, ['name']);

	return result;
},

isApp : function (appName, appGroup) {
// RETURN true IF AN APP EXISTS IN this.apps[appGroup]

	if ( this.apps[appGroup] == null )	return false;
	for ( var i in this.apps[appGroup] )
	{
		var app = this.apps[appGroup][i];
		if ( app.name.toLowerCase() == appName.toLowerCase() )
		{
			return true;
		}
	}

	return false;
},

getAppType : function (appName, appGroup) {
// RETURN THE TYPE OF AN APP GIVEN ITS NAME


	if ( this.apps[appGroup] == null )	return false;
	for ( var i in this.apps[appGroup] )
	{
		var app = this.apps[appGroup][i];
		if ( app.name.toLowerCase() == appName.toLowerCase() )
		{
			return app.type;
		}
	}

	return null;
},

// PARAMETER METHODS
addParameter : function (parameterObject) {
// ADD A PARAMETER OBJECT TO THE this.parameters ARRAY


	// REMOVE THE PARAMETER OBJECT IF IT EXISTS ALREADY
	var result = this._removeObjectFromArray(this.parameters, parameterObject, ['appname', 'name']);

	// DO THE ADD
	var requiredKeys = [ 'appname', 'name' ];
	return result = this._addObjectToArray(this.parameters, parameterObject, requiredKeys);
},


removeParameter : function (parameterObject) {
// REMOVE AN PARAMETER OBJECT FROM THE this.parameters ARRAY 


	var result = this._removeObjectFromArray(this.parameters, parameterObject, ['appname', 'name']);

	// RETURN TRUE OR FALSE
	return result;
},

isParameter : function (appName, parameterName) {
// RETURN true IF AN PARAMETER EXISTS IN this.parameters


	var parameters = getParametersByAppname(appName);
	if ( parameters == null )	return false;

	for ( var i in parameters )
	{
		if ( parameters[i].name.toLowerCase() == parameterName.toLowerCase() )
		{
			return true;
		}
	}

	return false;
},

getParametersByAppname : function (appname) {
// RETURN AN ARRAY OF PARAMETERS FOR THE GIVEN APPLICATION


	if ( appname == null )	return null;	
	var parameters = new Array;
	var params = dojo.clone(this.parameters);
	dojo.forEach(params, function(parameter)
	{
		if ( parameter.appname == appname ) parameters.push(parameter);
	});

	return parameters;
},

getParameter : function (appName, parameterName) {
// RETURN A NAMED PARAMETER FOR A GIVEN APPLICATION
// E.G., WHEN RETURNING VALUE TO DEFAULT


	if ( appName == null )	return null;
	if ( parameterName == null )	return null;

	var parameters = getParametersByAppname(appName);
	if ( parameters == null )	return false;

	for ( var i in parameters )
	{
		if ( parameters[i].name.toLowerCase() == parameterName.toLowerCase() )
		{
			return parameters[i];
		}
	}

	return null;
},


// STAGE METHODS
addStage : function (stageObject) {
// ADD AN APP OBJECT TO this.stages AND COPY ITS PARAMETERS
// INTO this.parameters.
// ALSO UPDATE THE REMOTE DATABASE.

	// DO THE ADD
	var result = this._addStage(stageObject);

	// ADD PARAMETERS FOR THIS STAGE TO this.stageparameters
	if ( result == true )	result = this.addStageParametersForStage(stageObject);

	if ( result == false )
	{
		return;
	}
	// ADD STAGE IN REMOTE DATABASE
	var url = Agua.cgiUrl + "workflow.cgi";
	var query = stageObject;
	query.username = Agua.cookie('username');
	query.sessionId = Agua.cookie('sessionId');
	query.mode = "addStage";

	this.doPut({ url: url, query: query });

	// RETURN TRUE OR FALSE
	return result;
},



_addStage : function (stageObject) {
// ADD A STAGE TO this.stages AND COPY ITS PARAMETERS
// INTO this.stageparameters 


	// ADD THE STAGE
	var requiredKeys = ['project', 'workflow', 'name', 'number'];
	var result = this._addObjectToArray(this.stages, stageObject, requiredKeys);
	if ( result == false )
	{
		return false;
	}

	//// ADD THE STAGE PARAMETERS
	//var stageParameters = this.getStageParameters(stageObject);
	//for ( var i = 0; i < stageParameters.length; i++ )
	//{
	//	this._addStageParameter(stageParameters[i]);
	//}

	return result;
},



insertStage : function (stageObject) {
// INSERT AN APP OBJECT INTO THE this.stages ARRAY,
// INCREMENTING THE number OF ALL DOWNSTREAM STAGES
// BEFOREHAND. DO THE SAME FOR THE this.stageparameters
// ENTRIES FOR THIS APP.
// THEN, MIRROR ON THE REMOTE DATABASE.


	//// INCREMENT STAGE NUMBERS OF DOWNSTREAM STAGES
	//this.incrementDownstreamStages(stageObject);

	// SANITY CHECK
	if ( stageObject == null )	return;

	// GET THE STAGES FOR THIS PROJECT AND WORKFLOW
	// ORDERED BY number
	var stages = this.getStagesByWorkflow(stageObject.project, stageObject.workflow);
	stages = this.sortHasharray(stages, "number");

	// GET THE INSERTION INDEX OF THE STAGE
	var index = stageObject.number - 1;

	// INCREMENT THE appnumber IN ALL DOWNSTREAM STAGES IN this.stageparameters
	var result = true;
	for ( var i = stages.length - 1; i > index - 1; i-- )
	{
		// GET THE STAGE PARAMETERS FOR THIS STAGE
		var stageParameters = this.getStageParameters(stages[i]);

		if ( stageParameters == null )
		{
			result = false;
		}

		// REMOVE EACH STAGE PARAMETER AND RE-ADD ITS UPDATED VERSION
		var thisObject = this;
		for ( var j = 0; j < stageParameters.length; j++ )
		{
			// REMOVE EXISTING STAGE
			if ( thisObject._removeStageParameter(stageParameters[j]) == false )
			{
				result = false;
			}

			// INCREMENT STAGE NUMBER
			stageParameters[j].appnumber = (i + 2).toString();

			// ADD BACK STAGE
			if ( thisObject._addStageParameter(stageParameters[j]) == false )
			{
				result = false;
			}				
		}

		//  ******** DEBUG ONLY ************
		//  ******** DEBUG ONLY ************
		//var updatedStageParameters = this.getStageParameters(stages[i]);
		//  ******** DEBUG ONLY ************
		//  ******** DEBUG ONLY ************
	}
	if ( result == false )
	{
		return false;
	}


	// INCREMENT THE number OF ALL DOWNSTREAM STAGES IN this.stages
	// NB: THE SERVER SIDE UPDATES ARE DONE AUTOMATICALLY
	for ( var i = stages.length - 1; i > index - 1; i-- )
	{
		// REMOVE FROM this.stages
		if ( this._removeStage(stages[i]) == false )
		{
			return false;
		}

		// INCREMENT STAGE NUMBER
		stages[i].number = (i + 2).toString();

		// ADD BACK TO this.stages
		if ( this._addStage(stages[i]) == false )
		{
			return false;
		}
	}

	////  ******** DEBUG ONLY ************
	////  ******** DEBUG ONLY ************
	//var newStages = this.getStagesByWorkflow(stageObject.project, stageObject.workflow);
	//newStages = this.sortHasharray(newStages, "number");
	//for ( var i = 0; i < newStages.length; i++ )
	//{
	//}
	////  ******** DEBUG ONLY ************
	////  ******** DEBUG ONLY ************


	// INSERT THE NEW STAGE (NO EXISTING STAGES WILL HAVE ITS number)
	// (NB: this._addStage CHECKS FOR REQUIRED FIELDS)
	result = this._addStage(stageObject);
	if ( result == false )
	{
		return false;
	}
	////  ******** DEBUG ONLY ************
	////  ******** DEBUG ONLY ************
	//for ( var i = 0; i < this.stages.length; i++ )
	//{
	//}
	////  ******** DEBUG ONLY ************
	////  ******** DEBUG ONLY ************




	// ADD PARAMETERS FOR THIS STAGE TO this.stageparameters
	if ( result == true ) result = this.addStageParametersForStage(stageObject);
	if ( result == false )
	{
		return false;
	}

	////  ******** DEBUG ONLY ************
	////  ******** DEBUG ONLY ************
	//for ( var i = 0; i < this.stageparameters.length; i++ )
	//{
	//}
	////  ******** DEBUG ONLY ************
	////  ******** DEBUG ONLY ************




	// ADD STAGE TO stage TABLE IN REMOTE DATABASE
	var url = Agua.cgiUrl + "workflow.cgi";
	var query = stageObject;
	query.username = Agua.cookie('username');
	query.sessionId = Agua.cookie('sessionId');
	query.mode = "insertStage";

	this.doPut({ url: url, query: query });

	// RETURN TRUE OR FALSE
	return true;
},


removeStage : function (stageObject) {
// REMOVE AN APP OBJECT FROM THE this.stages ARRAY
// AND SIMULTANEOUSLY REMOVE CORRESPONDING ENTRIES IN
// this.parameters FOR THIS STAGE.
// ALSO, MAKE removeStage CALL TO REMOTE DATABASE.


	// DO THE ADD
	var result = this._removeStage(stageObject);

	// ADD PARAMETERS FOR THIS STAGE TO this.stageparameters
	if( result == true ) 	result = this.removeStageParameters(stageObject);

	if ( result == false )
	{
		return;
	}

	// ADD STAGE TO stage TABLE IN REMOTE DATABASE
	var url = Agua.cgiUrl + "workflow.cgi";
	var query = stageObject;
	query.username = Agua.cookie('username');
	query.sessionId = Agua.cookie('sessionId');
	query.mode = "removeStage";

	this.doPut({ url: url, query: query });

	// RETURN TRUE OR FALSE
	return result;
},


_removeStage : function (stageObject) {
// REMOVE AN APP OBJECT FROM THE this.stages ARRAY
// AND SIMULTANEOUSLY REMOVE CORRESPONDING ENTRIES IN
// this.parameters FOR THIS STAGE

	var uniqueKeys = ['project', 'workflow', 'name', 'number'];
	var result = this._removeObjectFromArray(this.stages, stageObject, uniqueKeys);
	if ( result == false )
	{
		return false;
	}
	return result;
},


spliceStage : function (stageObject) {
// SPLICE OUT A STAGE FROM this.stages AND DECREMENT THE
// number OF ALL DOWNSTREAM STAGES.
// DO THE SAME FOR THE CORRESPONDING ENTRIES IN this.stageparameters


	// SANITY CHECK
	if ( stageObject == null )	return;

	// GET THE STAGES FOR THIS PROJECT AND WORKFLOW
	var stages = this.getStagesByWorkflow(stageObject.project, stageObject.workflow);

	// ORDER BY number AND SPLICE OUT THE STAGE 
	stages = this.sortHasharray(stages, "number");
	var index = stageObject.number - 1;
	stages.splice(index, 1);


	// REMOVE THE STAGE FROM THE this.stages
	var result = this._removeStage(stageObject);
	if ( result == false )
	{
		return false;
	}

	// REMOVE THE STAGE FROM this.stageparameters
	result = this.removeStageParameters(stageObject);
	if ( result == false )
	{
		return false;
	}

	// DECREMENT THE number OF ALL DOWNSTREAM STAGES IN this.stages
	// NB: THE SERVER SIDE UPDATES ARE DONE AUTOMATICALLY
	for ( var i = index; i < stages.length; i++ )
	{
		// REMOVE IT FROM THE this.stages ARRAY
		this._removeStage(stages[i]);
		stages[i].number = (i + 1).toString();
		this._addStage(stages[i]);

	}

	//  ******** DEBUG ONLY ************
	//  ******** DEBUG ONLY ************
	var newStages = this.getStagesByWorkflow(stageObject.project, stageObject.workflow);
	for ( var i = 0; i < newStages.length; i++ )
	{
	}

	//  ******** DEBUG ONLY ************
	//  ******** DEBUG ONLY ************


	// DECREMENT THE appnumber IN ALL DOWNSTREAM STAGES IN this.stageparameters



	////////  ******** DEBUG ONLY ************
	////////  ******** DEBUG ONLY ************
	//////for ( var i = 0; i < this.stageparameters.length; i++ )
	//////{
	//////}
	////////  ******** DEBUG ONLY ************
	////////  ******** DEBUG ONLY ************


	for ( var i = stages.length - 1; i > index - 1; i-- )
	{
		// REINCREMENT STAGE NUMBER TO GET ITS STAGE
		// PARAMETERS, WHICH HAVE NOT BEEN DECREMENTED YET
		stages[i].number = (i + 2).toString();


		// GET THE STAGE PARAMETERS FOR THIS STAGE
		var stageParameters = this.getStageParameters(stages[i]);
		if ( stageParameters == null )
		{
			return false;
		}

		// REMOVE EACH STAGE PARAMETER AND RE-ADD ITS UPDATED VERSION
		var thisObject = this;
		//dojo.forEach( stageParameters, function(stageParameter, i)
		for ( var j = 0; j < stageParameters.length; j++ )
		{
			// REMOVE EXISTING STAGE
			if ( thisObject._removeStageParameter(stageParameters[j]) == false )
			{
				result = false;
			}

			// DECREMENT STAGE NUMBER
			stageParameters[j].appnumber = (i + 1).toString();

			// ADD BACK STAGE
			if ( thisObject._addStageParameter(stageParameters[j]) == false )
			{
				result = false;
			}
		}

		// REDECREMENT STAGE NUMBER
		stages[i].number = (i + 1).toString();

	}

	////////  ******** DEBUG ONLY ************
	////////  ******** DEBUG ONLY ************
	//////
	//////for ( var i = 0; i < this.stageparameters.length; i++ )
	//////{
	//////}
	//////
	////////  ******** DEBUG ONLY ************
	////////  ******** DEBUG ONLY ************


	if ( result == false )
	{
		return false;
	}




	// REMOVE FROM REMOTE DATABASE DATABASE:
	var url = Agua.cgiUrl + "workflow.cgi";

	// GENERATE QUERY JSON FOR THIS WORKFLOW IN THIS PROJECT
	var query = stageObject;
	query.username = Agua.cookie('username');
	query.sessionId = Agua.cookie('sessionId');
	query.mode = "removeStage";

	this.doPut({ url: url, query: query});
},



isStage : function (stageName) {
// RETURN true IF AN APP EXISTS IN this.stages

	if ( this.stages == null )	return false;
	for ( var i in this.stages )
	{
		var stage = this.stages[i];
		if ( stage.name.toLowerCase() == stageName.toLowerCase() )
		{
			return true;
		}
	}

	return false;
},


getStageType : function (stageName) {
// RETURN true IF AN APP EXISTS IN this.stages

	if ( this.stages == null )	return false;
	for ( var i in this.stages )
	{
		var stage = this.stages[i];
		if ( stage.name.toLowerCase() == stageName.toLowerCase() )
		{
			return stage.type;
		}
	}

	return null;
},


getStages : function () {

	return dojo.clone(this.stages);
},


getStagesByWorkflow : function (project, workflow) {
// RETURN AN ARRAY OF STAGE HASHES FOR THIS PROJECT AND WORKFLOW

	if ( project == null )	return;
	if ( workflow == null )	return;

	var stages = dojo.clone(this.stages);
	//var temp = dojo.clone(stages);
	//stages = [];
	//stages = temp;

	var keyArray = ["project", "workflow"];
	var valueArray = [project, workflow];
	stages = this.filterByKeyValues(stages, keyArray, valueArray);


	return stages;
},



// STAGEPARAMETER METHODS
_addStageParameter : function (stageParameterObject) {
// ADD A STAGE PARAMETER OBJECT TO THE this.stageparameters ARRAY,
// REQUIRE THAT UNIQUE KEYS ARE DEFINED


	// DO THE ADD
	var uniqueKeys = ['username', 'project', 'workflow', 'appname', 'appnumber', 'name'];
	return this._addObjectToArray(this.stageparameters, stageParameterObject, uniqueKeys);
},


addStageParameter : function (stageParameterObject) {
// 1. REMOVE ANY EXISTING STAGE PARAMETER (UNIQUE KEYS: appname, appnumber, name)
// 		(I.E., THE SAME AS UPDATING AN EXISTING STAGE PARAMETER)
// 2. ADD A STAGE PARAMETER OBJECT TO THE this.stageparameters ARRAY
//   
// 3. ADD TO stageparameter TABLE IN REMOTE DATABASE


	// DO THE REMOVE
	this._removeStageParameter(stageParameterObject);

	// DO THE ADD
	var result = this._addStageParameter(stageParameterObject);
	if ( result == false )	return result;

	// REMOVE FROM REMOTE DATABASE DATABASE:
	// SET URL, ADD RANDOM NUMBER TO DISAMBIGUATE BETWEEN CALLS BY DIFFERENT
	// METHODS TO THE SERVER
	var url = this.cgiUrl + "workflow.cgi";

	// GENERATE QUERY JSON FOR THIS WORKFLOW IN THIS PROJECT
	var query = stageParameterObject;
	query.username = this.cookie('username');
	query.sessionId = this.cookie('sessionId');
	query.mode = "addStageParameter";

	this.doPut({ url: url, query: query});

	// RETURN TRUE OR FALSE
	return result;
},


_removeStageParameter : function (stageParameterObject) {
// REMOVE A stage PARAMETER OBJECT FROM THE this.stageparameters ARRAY,
// IDENTIFY OBJECT USING UNIQUE KEYS


	var uniqueKeys = ['username', 'project', 'workflow', 'appname', 'appnumber', 'name'];
	return this._removeObjectFromArray(this.stageparameters, stageParameterObject, uniqueKeys);
},


removeStageParameter : function (stageParameterObject) {
// 1. REMOVE A stage PARAMETER OBJECT FROM THE this.stageparameters ARRAY
// 2. REMOVE FROM stageparameter TABLE IN REMOTE DATABASE


	var result = this._removeStageParameter(stageParameterObject);
	if ( result == false )	return result;

	// REMOVE FROM REMOTE DATABASE:
	// SET URL, ADD RANDOM NUMBER TO DISAMBIGUATE BETWEEN CALLS BY DIFFERENT
	// METHODS TO THE SERVER
	var url = this.cgiUrl + "workflow.cgi";

	// GENERATE QUERY JSON FOR THIS WORKFLOW IN THIS PROJECT
	var query = stageParameterObject;
	query.username = this.cookie('username');
	query.sessionId = this.cookie('sessionId');
	query.mode = "addStageParameter";

	this.doPut({ url: url, query: query, timeout : 15000 });

	// RETURN TRUE OR FALSE
	return result;
},


getStageParameters : function (stageObject) {
// RETURN AN ARRAY OF STAGE PARAMETER HASHARRAYS FOR THE GIVEN STAGE


	var keys = ["project", "workflow", "name", "number"];
	var notDefined = this.notDefined (stageObject, keys);
	if ( notDefined.length != 0 )
	{
		return;
	}

	var stageParameters = dojo.clone(this.stageparameters);
	var keyArray = ["project", "workflow", "appname", "appnumber"];
	var valueArray = [stageObject.project, stageObject.workflow, stageObject.name, stageObject.number];
	stageParameters = this.filterByKeyValues(stageParameters, keyArray, valueArray);


	return stageParameters;
},

addStageParametersForStage : function (stageObject) {
// ADD this.parameters ENTRIES FOR A STAGE TO this.stageparameters


	if ( stageObject.name == null )	return null;
	if ( stageObject.number == null )	return null;

	// GET APP PARAMETERS	
	var parameters;
	if ( stageObject.owner == Agua.cookie('name') )
	{
		parameters = dojo.clone(this.getParametersByAppname(stageObject.name));
	}
	else {
		parameters = dojo.clone(this.getSharedParametersByAppname(stageObject.name, ));
	}


	// ADD STAGE project, workflow, AND number TO PARAMETERS
	dojo.forEach(parameters, function(parameter)
	{
		parameter.owner = stageObject.owner;
		parameter.project = stageObject.project;
		parameter.workflow = stageObject.workflow;
		parameter.appnumber = stageObject.number;
		parameter.appname = stageObject.name;
	});

	// ADD PARAMETERS TO this.stageparameters ARRAY
	var uniqueKeys = ['owner', 'project', 'workflow', 'appname', 'appnumber', 'name'];
	var addOk = true;
	var aguaObject = this;
	dojo.forEach(parameters, function(parameter)
	{
		if ( aguaObject._addObjectToArray(aguaObject.stageparameters, parameter, uniqueKeys) == false)
		{
			addOk = false;
		}
	});

	// ************** DEBUG ONLY ****************
	// ************** DEBUG ONLY ****************
	////for ( var i = 0; i < this.stageparameters.length; i++ )
	////{
	////}
	// ************** DEBUG ONLY ****************
	// ************** DEBUG ONLY ****************

	if ( ! addOk )
	{
		return;
	}

	return addOk;

},

removeStageParameters : function (stageObject) {
// REMOVE STAGE PARAMETERS FOR A STAGE FROM this.stageparameters 


	//////// ******************** DEBUG *********************
	//////// ******************** DEBUG *********************


	////////var header = '';
	////////for ( var name in this.stageparameters[0] )	{	header += name + "\t";	}
	////////
	////////for ( var i = 0; i < this.stageparameters.length; i++ )
	////////{
	////////	var info = '';
	////////	for ( var name in this.stageparameters[i] )
	////////	{
	////////		if ( this.stageparameters[i][name] != '' )
	////////			info += this.stageparameters[i][name] + "\t";
	////////		else info += "*\t";
	////////	}
	////////}


	//////// ******************** DEBUG *********************
	//////// ******************** DEBUG *********************

	if ( stageObject.name == null )	return null;
	if ( stageObject.number == null )	return null;

	// RETRIEVE PARAMETER NAMES FROM this.parameters
	var parameters = dojo.clone(this.getParametersByAppname(stageObject.name));

	// ADD STAGE project, workflow, AND number TO PARAMETERS
	var thisObject = this;
	dojo.forEach(parameters, function(parameter)
	{
		parameter.username = thisObject.cookie('username');
		parameter.project = stageObject.project;
		parameter.workflow = stageObject.workflow;
		parameter.appnumber = stageObject.number;
		parameter.appname = stageObject.name;
	});

	////// NB: !!! A SIMPLE FOR LOOP DOESN'T FIX THE VALUES IN THE PARAMETER !!!
	//////for ( var i = 0; i < parameters; i++ )
	//////{
	//////	parameters[i].username = this.cookie('username');
	//////	parameters[i].project = stageObject.project;
	//////	parameters[i].workflow = stageObject.workflow;
	//////	parameters[i].appnumber = stageObject.number;
	//////	parameters[i].appname = stageObject.name;
	//////}


	// REMOVE PARAMETERS FROM this.stageparameters
	var uniqueKeys = [ "username", "project", "workflow", "appname", "appnumber", "name"];
	var removeOk = true;
	for ( var j = 0; j < parameters.length; j++ )
	{

		var notDefined = this.notDefined (parameters[j], uniqueKeys);
		if ( notDefined.length > 0 )
		{
			removeOk = false;
			continue;
		}

		var removeSuccess = this._removeObjectFromArray(this.stageparameters, parameters[j], uniqueKeys);

		////////// ***** DEBUG ONLY *****
		////////if ( removeSuccess == false )
		////////{				
		////////	removeOk = false;
		////////}
	}

	if ( removeOk == false )
	{
		return false;
	}

	return true;
},

getStageParametersByApp : function (appname) {
// RETURN AN ARRAY OF PARAMETER HASHARRAYS FOR THE GIVEN APPLICATION


	var stageParameters = new Array;
	dojo.forEach(this.stageparameters, function(stageparameter)
	{
		if ( stageparameter.appname == appname )	stageParameters.push(stageparameter);
	});


	return stageParameters;
},

getParameterValidity : function (stageParameterObject, booleanValue) {
// GET THE BOOLEAN parameterValidity VALUE FOR A STAGE PARAMETER


	if ( booleanValue != null )
	{
		return null;
	}

	var isValid = this._parameterValidity(stageParameterObject, booleanValue);

	return isValid;
},

setParameterValidity : function (stageParameterObject, booleanValue) {
// SET THE BOOLEAN parameterValidity VALUE FOR A STAGE PARAMETER


	if ( booleanValue == null )
	{
		return null;
	}

	return this._parameterValidity(stageParameterObject, booleanValue);
},

_parameterValidity : function (stageParameterObject, booleanValue) {
// RETURN THE parameterValidity BOOLEAN FOR A STAGE PARAMETER
// OR SET IT IF A VALUE IS SUPPLIED


	//var filtered = dojo.clone(this.stageparameters);
	//var keys = ["appname"];
	//var values = ["image2eland.pl"];
	//filtered = this.filterByKeyValues(filtered, keys, values);

	var uniqueKeys = ['username', 'project', 'workflow', 'appname', 'appnumber', 'name'];
	var valueArray = new Array;
	for ( var i = 0; i < uniqueKeys.length; i++ )
	{
		valueArray.push(stageParameterObject[uniqueKeys[i]]);
	}
	var stageParameter = this.getEntry(this.stageparameters, uniqueKeys, valueArray);
	if ( stageParameter == null )
	{
		return null;
	}

	if ( booleanValue != null )
	{
		// SET isValid BOOLEAN VALUE
		stageParameter.isValid = booleanValue;		
		var success = this._removeStageParameter(stageParameter);
		if ( success == false )
		{
			return null;
		}


		success = this._addStageParameter(stageParameter);			
		if ( success == false )
		{
			return null;
		}

		return true;
	}
	else
	{
		return stageParameter.isValid;
	}
},


// REPORT METHODS
getReports : function () {
// RETURN A COPY OF THE this.reports ARRAY


	var stages = this.getStages();
	if ( stages == null )	return;

	var keys = [ "type" ];
	var values = [ "report" ];
	var reports = this.filterByKeyValues(stages, keys, values);

	return reports;
},

getReportsByProjectWorkflow : function (projectName, workflowName) {
// RETURN AN ARRAY OF REPORT HASHES FOR THE SPECIFIED PROJECT AND WORKFLOW

	var reports = this.getReports();
	if ( reports == null )	return null;

	var keyArray = ["project", "workflow"];
	var valueArray = [projectName, workflowName];
	return this.filterByKeyValues(reports, keyArray, valueArray);
},

removeReport : function (reportObject) {
// REMOVE A REPORT OBJECT FROM THE this.reports ARRAY


	// REMOVE REPORT FROM this.reports
	var requiredKeys = ["project", "workflow", "name"];
	var success = this._removeObjectFromArray(this.reports, reportObject, requiredKeys);

	// REMOVE REPORT FROM this.groupmembers IF PRESENT
	var groupNames = this.getGroupsByReport(reportObject.name);
	for ( var i = 0; i < groupNames.length; i++ )
	{
		if ( this.removeReportFromGroup(groupNames[i], reportObject) == false )
			success = false;
	}

	return success;
},

isReport : function (reportName) {
// RETURN true IF A REPORT EXISTS IN this.reports


	if ( this.reports == null )	return false;
	for ( var i in this.reports )
	{
		var report = this.reports[i];
		if ( report.name.toLowerCase() == reportName.toLowerCase() )
		{
			return true;
		}
	}

	return false;
},

addReport : function (reportObject) {
// ADD A REPORT 


	// DO THE ADD
	var requiredKeys = ["project", "workflow", "name"];
	var success = this._addObjectToArray(this.reports, reportObject, requiredKeys);

	return success;
},

// VIEW METHODS
getViewObject : function (projectName, viewName) {
// RETURN AN ARRAY OF VIEW HASHES FOR THE SPECIFIED PROJECT AND WORKFLOW
	var views = this.getViews();
	if ( views == null )	return [];
	var keyArray = ["project", "view"];
	var valueArray = [projectName, viewName];
	var views = this.filterByKeyValues(views, keyArray, valueArray);

	return views[0];
},

getViews : function () {
// RETURN A COPY OF THE this.views ARRAY

	if ( this.views == null )	return;

	var views = dojo.clone(this.views);

	return views;
},


removeView : function (viewObject) {
// REMOVE A VIEW OBJECT FROM THE this.views ARRAY


	// REMOVE VIEW FROM this.views
	var requiredKeys = ["project", "view"];
	var success = this._removeObjectFromArray(this.views, viewObject, requiredKeys);
	if ( success == true ) {
		return true;
	}
	else {
		return false;
	}
},

isView : function (projectName, viewName) {
// RETURN true IF A VIEW EXISTS FOR THE PARTICULAR PROJECT AND WORKFLOW


	var viewObjects = this.getViewsByProject(projectName);
	for ( var i in viewObjects )
	{
		var viewObject = viewObjects[i];
		if ( viewObject.view.toLowerCase() == viewName.toLowerCase() )
		{
			return true;
		}
	}

	return false;
},

addView : function (viewObject) {
// ADD A VIEW TO this.views AND SAVE ON REMOTE SERVER


	// DO THE ADD
	var requiredKeys = ["project", "view"];
	return this._addObjectToArray(this.views, viewObject, requiredKeys);
},

viewNames : function (projectName) {
// RETURN AN ARRAY OF ALL VIEW NAMES IN this.views

	var views = this.getViewsByProject(projectName);

	return this.hasharrayKeyToArray(views, "view");
},

getViewsByProject : function (projectName) {
// RETURN AN ARRAY OF VIEW HASHES FOR THE SPECIFIED PROJECT AND WORKFLOW

	var views = this.getViews();
	if ( views == null )	return [];

	var keyArray = ["project"];
	var valueArray = [projectName];
	return this.filterByKeyValues(views, keyArray, valueArray);
},



getViewSpecies : function (projectName, viewName) {
// GET THE UNIQUE SPECIES (AND BUILD) FOR A GIVEN VIEW


	if ( projectName == null || ! projectName )
	{
		return;
	}

	var viewfeatures = this.getViewFeatures(projectName, viewName);
	if ( viewfeatures == null || viewfeatures.length == 0 )	return new Array;

	var speciesHash = new Object;
	speciesHash.species = viewfeatures[0].species;
	speciesHash.build = viewfeatures[0].build;

	return speciesHash;
},


getViewFeatures : function (projectName, viewName) {
// GET THE UNIQUE SPECIES (AND BUILD) FOR A GIVEN VIEW


	if ( projectName == null || ! projectName )
	{
		return;
	}

	if ( this.viewfeatures == null )	return new Array;
	var viewfeatures = dojo.clone(this.viewfeatures);
	var keyArray = ["project", "view"];
	var valueArray = [projectName, viewName];
	viewfeatures = this.filterByKeyValues(viewfeatures, keyArray, valueArray);

	return viewfeatures;
},

getViewProjects : function () {

	if ( this.viewfeatures == null )	return [];
	var viewfeatures = dojo.clone(this.viewfeatures);
	var projects = this.hasharrayKeyToArray(viewfeatures, "project");
	projects = this.uniqueValues(projects);

	return projects;
},

// FEATURE METHODS
getViewProjectWorkflows : function (projectName) {
	if ( projectName == null || ! projectName )
	{
		return;
	}

	if ( this.features == null )	return new Array;
	var features = dojo.clone(this.features);
	var keyArray = ["project"];
	var valueArray = [projectName];
	features = this.filterByKeyValues(features, keyArray, valueArray);
	var workflows = new Array;
	for ( var i = 0; i < features.length; i++ )
		workflows.push(features[i].workflow);

	workflows = this.uniqueValues(workflows);

	return workflows;
},

getViewWorkflowFeatures : function (projectName, workflowName) {

	if ( this.features == null || this.features.length == 0 )	return new Array;
	var features = dojo.clone(this.features);
	var keyArray = ["project", "workflow"];
	var valueArray = [projectName, workflowName];
	features = this.filterByKeyValues(features, keyArray, valueArray);

	return features;
},

getViewSpeciesFeatureNames : function (projectName, workflowName, speciesName, buildName) {
	// GET THE FEATURE NAMES FOR A GIVEN PROJECT, WORKFLOW AND SPECIES BUILD

	var features = this.getViewSpeciesFeatures(projectName, workflowName, speciesName, buildName);
	if ( this.features == null || this.features.length == 0 )	return new Array;
	var featureNames = new Array;
	for ( var i = 0; i < features.length; i++ )
		featureNames.push(features[i].feature);

	featureNames = this.uniqueValues(featureNames);

	return featureNames;
},

getViewSpeciesFeatures : function (projectName, workflowName, speciesName, buildName) {
	// GET THE FEATURES FOR A GIVEN PROJECT, WORKFLOW AND SPECIES BUILD

	if ( projectName == null || ! projectName
		|| workflowName == null || ! workflowName
		|| speciesName == null || ! speciesName )
	{
		return;
	}

	if ( this.features == null || this.features.length == 0 )	return new Array;
	var features = dojo.clone(this.features);
	var keyArray = ["project", "workflow", "species", "build"];
	var valueArray = [projectName, workflowName, speciesName, buildName];
	features = this.filterByKeyValues(features, keyArray, valueArray);

	return features;
},








removeViewFeature : function (featureObject) {

	// REMOVE LOCALLY AND THEN ON THE REMOTE
	var requiredKeys = ["project", "view", "feature"];
	var success = this._removeObjectFromArray(this.viewfeatures, featureObject, requiredKeys);
	if ( success )
	else
	{
		return;
	}
	var url = Agua.cgiUrl + "workflow.cgi";
	featureObject.username = Agua.cookie('username');
	featureObject.sessionId = Agua.cookie('sessionId');
	featureObject.mode = "removeViewFeature";

	this.doPut({ url: url, query: featureObject });
},

isFeature : function (featureObject) {
// RETURN true IF THE FEATURE ALREADY EXISTS IN THE VIEW


	var viewfeatures = this.getViewFeatures(featureObject.project, featureObject.view);
	if ( this._objectInArray(viewfeatures, featureObject, ["project", "view", "feature", "species", "build"]))
	{
		return true;
	}

	return false;
},

addViewFeature : function (featureObject) {

	// CHECK IF FEATURE ALREADY EXISTS
	if ( this.isFeature(featureObject) )
		return;

	// ADD LOCALLY 
	var requiredKeys = ["project", "view", "feature", "species", "build"];
	return this._addObjectToArray(this.viewfeatures, featureObject, requiredKeys);
},



// CLUSTER METHODS
getClusterObject : function (projectName, clusterName) {
// RETURN AN ARRAY OF CLUSTER HASHES FOR THE SPECIFIED PROJECT AND WORKFLOW
	var clusters = this.getClusters();
	if ( clusters == null )	return [];
	var keyArray = ["project", "cluster"];
	var valueArray = [projectName, clusterName];
	var clusters = this.filterByKeyValues(clusters, keyArray, valueArray);

	return clusters[0];
},

getClusters : function () {
// RETURN A COPY OF THE this.clusters ARRAY

	if ( this.clusters == null )	return;

	var clusters = dojo.clone(this.clusters);

	return clusters;
},

removeCluster : function (clusterObject) {
// REMOVE A CLUSTER OBJECT FROM THE this.clusters ARRAY


	// REMOVE CLUSTER FROM this.clusters
	var requiredKeys = ["cluster", "instancetype", "minnodes", "maxnodes"];
	var success = this._removeObjectFromArray(this.clusters, clusterObject, requiredKeys);
	if ( success == true ) {
		return true;
	}
	else {
		return false;
	}
},

isCluster : function (clusterName) {
// RETURN true IF A CLUSTER EXISTS


	var clusterObjects = this.getClusters();
	var inArray = this._objectInArray(clusterObjects, { cluster: clusterName }, ["cluster"]);	

	return inArray;
},

addCluster : function (clusterObject) {
// ADD A CLUSTER TO this.clusters AND SAVE ON REMOTE SERVER


	// DO THE ADD
	var requiredKeys = ["username", "cluster", "minnodes", "maxnodes", "instancetype"];
	return this._addObjectToArray(this.clusters, clusterObject, requiredKeys);
},


// HOUSEKEEPING METHODS
hideLoader : function() {
/* HIDE LOADING SCREEN


	dojo.fadeOut({

		node:"splashNode",
		duration:300,
		onEnd: function(){
			dojo.style("splashNode", "display", "none");
		}
	}).play();

	dojo.fadeOut({

		node:"backgroundNode",
		duration:300,
		onEnd: function(){
			dojo.style("backgroundNode", "display", "none");
		}
	}).play();
*/
},

showPreloadMessage : function (message) {
// SHOW LOADING PROGRESS

//		if ( message == null || ! message )
//		{
//			message = '';
//		}
//		
//return;
//
//		
//		this.messageNode.innerHTML += message + "<br>";
},

cookie : function (name, value) {
// SET OR GET COOKIE-CONTAINED USER ID AND SESSION ID


	if ( value != null )
	{
		this.cookies[name] = value;
	}
	else if ( name != null )
	{
		return this.cookies[name];
	}


	return 0;
},

loadCSSFile : function (cssFile) {
// LOAD A CSS FILE IF NOT ALREADY LOADED, REGISTER IN this.loadedCssFiles


	if ( this.loadedCssFiles == null || ! this.loadedCssFiles )
	{
		this.loadedCssFiles = new Object;
	}

	if ( ! this.loadedCssFiles[cssFile] )
	{

		var cssNode = document.createElement('link');
		cssNode.type = 'text/css';
		cssNode.rel = 'stylesheet';
		cssNode.href = cssFile;
		document.getElementsByTagName("head")[0].appendChild(cssNode);

		this.loadedCssFiles[cssFile] = 1;
	}
	else
	{
	}


	return this.loadedCssFiles;
}





}); // end of Agua

dojo.addOnLoad(
	function()
	{
	}
);

