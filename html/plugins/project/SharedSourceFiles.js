
dojo.provide("plugins.project.SharedSourceFiles");

// DISPLAY THE USER'S OWN SOURCES DIRECTORY AND ALLOW
// THE USER TO BROWSE AND MANIPULATE FOLDERS AND FILES

// LATER FOR MENU: HOW TO DYNAMICALLY
// ENABLE / DISABLE MENU ITEM
// attr('disabled', bool) 


// INHERITS
dojo.require("plugins.project.ProjectFiles");


dojo.declare( "plugins.project.SharedSourceFiles",
	[ plugins.project.ProjectFiles ],
{
	//Path to the template of this widget. 
	templatePath: dojo.moduleUrl("plugins", "project/templates/filesystem.html"),

	// Calls dijit._Templated.widgetsInTemplate
	widgetsInTemplate : true,

	// PARENT WIDGET
	parentWidget : null,

	// PROJECT NAME AND WORKFLOW NAME IF AVAILABLE
	project : null,

	// DEFAULT TIME (milliseconds) TO SLEEP BETWEEN FILESYSTEM LOADS
	sleep : 250,

	// CONSTRUCTOR	
	constructor : function(args)
	{		

		// SET PANE ID
		this.paneId = args.paneId;

		// SET ARGS
		this.project = args.project;
		this.parentWidget = args.parentWidget;
		this.attachNode = args.attachNode;

		// LOAD SORIA AND FILEPICKER CSS
        this.loadCSS();
	},

	// DO STARTUP
	postCreate: function()
	{

		this.startup();
	},

	// START MENUS
	startup : function ()
	{

		// SET UP THE ELEMENT OBJECTS AND THEIR VALUE FUNCTIONS
		this.inherited(arguments);

		//// ADD THE PANE TO THE TAB CONTAINER
		//this.attachNode.appendChild(this.mainNode);
		//
		//// LOAD APPLICATIONS DND SOURCE AND TARGET INTO FRAMEWORK
		//this.load();
	},

	// RELOAD THE WHOLE PROJECT TAB
	reload : function ()
	{
		this.load();
	},

	// GET DIRECTORIES TO SEARCH FOR FILES
	getDirectories : function ()
	{
		return Agua.getSharedSources();
	},


	setFileDrag : function (directory)
	{

		var owner = directory.owner;
		var name = directory.name;
		var description = directory.description;
		if ( ! description ) { description = '' };

		var titlePane = this.createTitlePane(
		{
			owner: owner,
			name: name,
			description: description,
			open: this.open
		});

		var thisStore = this.createStore(directory);

		// GENERATE NEW FileDrag OBJECT
		var fileDrag = new plugins.files.FileDrag(
			{
				style: "height: auto; width: 100%; minHeight: 50px;",
				store: thisStore,
				fileMenu: this.fileMenu,
				folderMenu: this.folderMenu,
				workflowMenu: this.workflowMenu,
				owner: owner,
				parentWidget: this
			}
		);

		// SET PATH AS LOCATION
		fileDrag.path = directory.location;                    

		// START UP FileDrag
		fileDrag.startup();

		// ADD fileDrag TO TITLE PANE
		titlePane.containerNode.appendChild(fileDrag.domNode);

	},

	// CREATE URL FOR STORE
	createUrl : function (directory)
	{
		// SET URL
		var url = Agua.cgiUrl;
		url += "project.cgi?";
		url += "mode=fileSystem";
		url += "&sessionId=" + Agua.cookie('sessionId');
		url += "&requestor=" + Agua.cookie('username');
		url += "&username=" + directory.owner;
		url += "&groupname=" + directory.groupname;
		url += "&location=" + directory.location;

		return url;		
	},

	// CREATE TITLE PANE
	createTitlePane : function (args)
	{

		var titlePaneNode = document.createElement('div');
		this.rowsNode.appendChild(titlePaneNode);

		// CREATE TITLE TABLE
		var title = "<table class='titlePane'><tr><td class='owner'>" + args.owner + "</td><td class='name'><span class='type'>" + this.type + "</span><span class='name'>" + args.name + "</span></td><td class='description'>" + args.description + "</td></tr></table>";

		// CREATE TITLE PANE			
		var titlePane = new dijit.TitlePane(
			{
				title: title,
				open: this.open
			},
			titlePaneNode
		);

		return titlePane;
	},



	//// CREATE URL FOR STORE
	//createUrl : function (directory)
	//{
	//	// SET URL
	//	var url = Agua.cgiUrl;
	//	url += "project.cgi?";
	//	
	//	var query = "mode=fileSystem";
	//	query += "&username=" + Agua.cookie('username');
	//	query += "&sessionId=" + Agua.cookie('sessionId');
	//
	//	//// USE STRAIGHT LOCATION IF BEGINS WITH '/'
	//	//if ( directory.location.match(/^\//) )
	//	//{
	//		query += "&location=" + directory.location;
	//	//}
	//	//else
	//	//{
	//	//	// SET LOCATION TO BE USED INSTEAD OF FILEROOT
	//	//	// ASSUMES ALL SOURCES ARE AT LEAST TWO FOLDERS DEEP IN THE FILESYSTEM
	//	//	var sourceBase = directory.location.match(/^(.+)\/[^\/]+$/)[1];
	//	//	var sourceFolder = directory.location.match(/^.+\/([^\/]+)$/)[1];
	//	//	query += "&location=" + sourceBase;
	//	//}
	//
	//
	//	return url + query;
	//},


	// CREATE STORE FOR FILE DRAG
	createStore : function (directory)
	{


		// SET URL
		var url = this.createUrl(directory);

		// CREATE STORE
		var thisStore = new dojox.data.FileStore(
			{
				//id: paneNodeId + "-fileStore",
				url: url,
				pathAsQueryParam: true
			}
		);

		// SET FILE STORE path TO project
		thisStore.preamble = function()
		{
			this.store.path = this.arguments[0].path;                        
		};

		return thisStore;		
	}

});
