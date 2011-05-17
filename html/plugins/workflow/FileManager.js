dojo.provide("plugins.workflow.FileManager");

dojo.require("dojox.widget.RollingList");

// DISPLAY THE USER'S PROJECTS DIRECTORY AND ALLOW
// THE USER TO BROWSE FILES AND MANIPULATE WORKFLOW
// FOLDERS AND FILES

	// THE callback FUNCTION IS PASSED THROUGH FileManager, FileSelector
	// AND _GroupSelectorPane:
	//
	//
	//	FileManager:
	//		preamble: function()
	//		constructor : function(args)  ??? not necessary
	// 		loadProjectTab: function ()
	//
	//
	//	FileSelector
	//		constructor(args)		
	//		getPaneForItem: function(/*item*/ item, /* dijit._Contained */ parentPane, /* item[]? */ children){
	//f
	//
	//	_GroupSelectorPane
	//		constructor(args)		
	//		createMenu : function (type) ... dojo.connect(mItem1, "onClick", function()
	//
	//


dojo.require("plugins.files.FileDrag");

// FILESTORE
dojo.require("dojox.data.FileStore");

// DIALOGUE
dojo.require("dojox.widget.Dialog");

// DnD
dojo.require("dojo.dnd.Source"); // Source & Target
dojo.require("dojo.dnd.Moveable");
dojo.require("dojo.dnd.Mover");
dojo.require("dojo.dnd.move");

// TITLE PANE
dojo.require("dijit.TitlePane");

// TOOLTIP
dojo.require("dijit.Tooltip");

// MENUS
dojo.require("plugins.files.FileSelectorMenu");
dojo.require("plugins.files.WorkflowSelectorMenu");


// WidgetFramework ancestor
dojo.require("plugins.core.WidgetFramework");

dojo.declare( "plugins.workflow.FileManager",
	[ dojox.widget.RollingList, plugins.core.Common ],
{
// PANE IDENTIFIER
paneId : '',

// STORE DND SOURCE ID AND DND TARGET ID
sourceId: '',

// PROJECT AND WORKFLOW
projects : new Array,

// DIALOGUE TO DISPLAY FILE MANAGER
fileManagerDialog : null,

// OWNER workflowObject
workflowObject : null,

// ID FOR THIS FILE MANAGER DIALOG PANE
dialogId : null,

// CSS FILES
cssFiles : [ dojo.moduleUrl("plugins") + "/workflow/css/filemanager.css"],

// callback FUNCTION AND DATA FROM OBJECT THAT GENERATED THE FileManager
callback : null,

// CORE WORKFLOW OBJECTS
core : null,

	////}

preamble: function(){
	this.callback = arguments[0].callback;
},

// CONSTRUCTOR	
constructor : function(args) {

	this.core = args.core;

	this.selectCallback = args.selectCallback;
	this.addCallback = args.addCallback;

	if ( ! args.paneId )
	{
		return;
	}

	// GET PROJECT AND WORKFLOW FROM ARGS
	this.parentWidget = args.parentWidget;

	// SET Workflow.fileManager AS THIS WIDGET
	this.parentWidget.fileManager = this;

	// SET PANE ID
	this.paneId = args.paneId;

	// LOAD SORIA AND FILEPICKER CSS
	this.loadCSS();

	// SET MENUS
	this.setMenus();

	// INSTANTIATE this.fileManagerDialog AND LOAD ITS PANES
	this.loadProjectTab();

	// SET ORIGINATOR, I.E. WORKFLOW
	this.workflowObject = args.workflowObject;

	// SUBSCRIBE TO UPDATES
	Agua.updater.subscribe(this, "updateProjects");

	//// LOAD APPLICATIONS DND SOURCE AND TARGET INTO FRAMEWORK
	//this.loadProjectTab();

	//this.hide();
},


updateProjects : function (args) {
// RELOAD THE GROUP COMBO AND DRAG SOURCE AFTER CHANGES
// TO SOURCES OR GROUPS DATA IN OTHER TABS


	// SET DRAG SOURCE
	this.loadProjectTab();
},


setMenus : function () {

	this.fileMenu = new plugins.files.FileSelectorMenu(
		{
			type: "file",
			selectCallback: this.selectCallback,
			addCallback: this.addCallback
		}
	);
	this.folderMenu = new plugins.files.FileSelectorMenu(
		{
			type: "directory",
			selectCallback: this.selectCallback,
			addCallback: this.addCallback
		}
	);

	this.workflowMenu = new plugins.files.FileSelectorMenu(
		{
			type: "directory",
			selectCallback: this.selectCallback,
			addCallback: this.addCallback
		}
	);


	// STOP PROPAGATION TO NORMAL RIGHTCLICK CONTEXT MENU
	dojo.connect(this.fileMenu.menu.domNode, "oncontextmenu", function (event)
	{
		event.stopPropagation();
	});

	//// STOP PROPAGATION TO NORMAL RIGHTCLICK CONTEXT MENU
	//dojo.connect(this.folderMenu.menu.domNode, "oncontextmenu", function (event)
	//{
	//	event.stopPropagation();
	//});
	//
	//// STOP PROPAGATION TO NORMAL RIGHTCLICK CONTEXT MENU
	//dojo.connect(this.workflowMenu.menu.domNode, "oncontextmenu", function (event)
	//{
	//	event.stopPropagation();
	//});
},

hide : function () {

	if ( this.fileManagerDialog == null || !this.fileManagerDialog )
	{
		return;
	}

	// SET fileManagerDialog
	this.fileManagerDialog.domNode.setAttribute('style', 'visibility: hidden');

	dojo.style(this.fileManagerDialog.containerNode, {
		opacity: 0,
		overflow: "hidden"
	});
},


show : function (parameterWidget) {
	console.dir(parameterWidget.args);

	// SET fileManagerDialog
	this.fileManagerDialog.domNode.setAttribute('style', 'visibility: visible');

	dojo.style(this.fileManagerDialog.containerNode, {
		opacity: 1,
		overflow: "auto"
	});

	// SET parameterWidget IN THIS AND ITS MENUS
	this.parameterWidget = parameterWidget;
	this.fileMenu.parameterWidget = parameterWidget;
	this.folderMenu.parameterWidget = parameterWidget;
	this.workflowMenu.parameterWidget = parameterWidget;

	// SHOW THE DIALOG
	this.fileManagerDialog.show();

	// OPEN THE TITLE PANE FOR THIS PROJECT
	var project = this.parameterWidget.project;
	var workflow = this.parameterWidget.workflow;
	for ( var i = 0; i < this.projects.length; i++ )
	{
		var titlePane = this.projects[i];
		if ( titlePane.title == project )
		{
			// OPEN PROJECT PANE
			titlePane.open = false;
			titlePane.toggle();


			// SELECT WORKFLOW IF DEFINED
			if ( workflow != null && workflow )
				titlePane.fileSelector.selectChild(workflow);
		}
	}
},

disableMenus : function () {

	this.fileMenu.disable();
	this.folderMenu.disable();
	this.workflowMenu.disable();
},

enableMenus : function () {

	this.fileMenu.enable();
	this.folderMenu.enable();
	this.workflowMenu.enable();
},

getProjects : function() {
	var projects = Agua.getProjects();

	// COPY name TO project FIELD
	for ( var i = 0; i < projects.length; i++ )
	{
		projects[i].project = projects[i].name;
	}

	return projects;
},


getSources : function() {
	var sources = Agua.getSources();

	// COPY name TO source FIELD
	for ( var i = 0; i < sources.length; i++ )
	{
		sources[i].source = sources[i].name;
	}

	return sources;
},	



loadProjectTab : function () {

	var dialogId = dojo.dnd.getUniqueId();		
	this.dialogId = dialogId;

	// CREATE DIALOGUE
	var fileManagerDialog = new dojox.widget.Dialog ({
		id: dialogId,
		dimensions: [800,"auto"],
		draggable: true,
		align: 'left',
		sizeDuration: 200,
		sizeMethod:       "combine",
		viewportPadding : "125",
		showTitle: "true",
		title: "File Selector"

	}, "fileManagerDialog" );		

	fileManagerDialog.set('widgetId', dialogId + "_widget");



	// OVERRIDE DIALOG _position METHOD:
	// 1. MAKE SURE FILE WINDOW STAYS PUT ON PAGE SCROLL
	// 2. FIX DIALOG DIMENSIONS
	fileManagerDialog._position = function() {

		if(this._sizing){
			this._sizing.stop();	
			this.disconnect(this._sizingConnect);
		}


		// EXTRACTED FROM dijit.Dialog			
		// summary: Position modal dialog in the viewport. If no relative offset
		//		in the viewport has been determined (by dragging, for instance),
		//		center the node. Otherwise, use the Dialog's stored relative offset,
		//		and position the node to top: left: values based on the viewport.
		if(!dojo.hasClass(dojo.body(),"dojoMove")){

			var node = this.domNode;
			var viewport = dijit.getViewport();
			var p = this._relativePosition;
			var mb = p ? null : dojo.marginBox(node);
			dojo.style(node,{
				left: Math.floor(viewport.l + (p ? p.l : (viewport.w - mb.w) / 2)) + "px",
				top: "100px"
			});

		}		
		//this.inherited("_position", arguments);

		if(!this.open){ dojo.style(this.containerNode, "opacity", 0); }
		var pad = this.viewportPadding * 2; 

		var props = {
			node: this.domNode,
			duration: this.sizeDuration || dijit._defaultDuration,
			easing: this.easing,
			method: this.sizeMethod
		};

		// CHANGE _displaysize.h FROM "auto" TO 800, IGNORE this._vp 
		this._displaysize = { h: 0, w: 800 };
		var ds = this._displaysize;
		props['width'] = ds.w;
		props['height'] = ds.h;

		//var ds = this._displaysize;
		//props['width'] = ds.w = (ds.w + pad >= this._vp.w || this.fixedSize) 
		//	? this._vp.w - pad : ds.w;
		//	
		//props['height'] = ds.h = (ds.h + pad >= this._vp.h || this.fixedSize) 
		//	? this._vp.h - pad : ds.h;

		this._sizing = dojox.fx.sizeTo(props);
		this._sizingConnect = this.connect(this._sizing,"onEnd","_showContent");
		this._sizing.play();

	};


	// SET CLASS FOR STYLE OF INSERTED PANE
	dojo.attr(dojo.byId(dialogId), 'class', 'fileManager dijitDialog');

	// SET TITLE CLASS
	dojo.attr(fileManagerDialog.titleBar, 'class', 'title');

	var projects = this.getProjects();
	if ( projects == null )
	{
		return;
	}

	////////fileManagerDialog.setContent(paneNode);
	var paneNode = document.createElement('div');

	// ADD fileManager CLASS TO NODE
	dojo.addClass(paneNode, 'fileManager');

	// *****************
	// *** PROJECTS ****
	// *****************

	var projectContainerNode = document.createElement('div');
	paneNode.appendChild(projectContainerNode);
	for ( var i = 0; i < projects.length; i++ )
	{

		var project = projects[i].project;
		var owner = projects[i].owner;

		// CREATE TITLE 
		var titlePaneNode = document.createElement('div');
		paneNode.appendChild(titlePaneNode);
		//this.leftPane.domNode.appendChild(titlePaneNode);

		// CREATE TITLE PANE			
		var titlePane = new dijit.TitlePane(
			{
				title: project,
				open: false
			},
			titlePaneNode
		);

		var titleNodeId = dojo.dnd.getUniqueId();
		titlePane.titleNode.setAttribute( 'id', titleNodeId );

		var tooltipLabel = projects[i].description;
		var tooltip = new dijit.Tooltip(
			{
				label: tooltipLabel,
				connectId:[titleNodeId],
				position : ['above']
			}
		);

		// SET URL
		var url = Agua.cgiUrl;
		url += "project.cgi?";
		url += "mode=fileSystem";
		url += "&username=" + Agua.cookie('username');
		url += "&sessionId=" + Agua.cookie('sessionId');

		// CREATE STORE
		var selectorStore = new dojox.data.FileStore(
			{
				//id: paneNodeId + "-fileStore",
				url: url,
				pathAsQueryParam: true
			}
		);

		// SET FILE STORE path TO project
		selectorStore.preamble = function()
		{
			this.store.path = this.arguments[0].path;                        
		};

		// GENERATE NEW FileDrag OBJECT
		var fileSelector = new plugins.files.FileDrag(
			{
				style: "height: auto; width: 100%; minHeight: 50px;",
				store: selectorStore,
				fileMenu: this.fileMenu,
				folderMenu: this.folderMenu,
				workflowMenu: this.workflowMenu,
				owner: owner,
				parentWidget: this
			}
		);

		// SET PATH FOR THIS PROJECT
		fileSelector.path = project;                    		

		// START UP FileDrag
		fileSelector.startup();

		// SHOW CLOSE BUTTON
		fileManagerDialog.closeButtonNode.setAttribute('style', 'opacity : 1');

		// ADD fileSelector TO TITLE PANE
		titlePane.containerNode.appendChild(fileSelector.domNode);

		// SET TITLE PANE fileSelector
		titlePane.fileSelector = fileSelector;

		//// APPEND FileDrag TO PANE
		//dojo.byId(sourceId).appendChild(paneNode);

		this.projects.push(titlePane);

	} // for loop on projects.length		

/*
	//// AVOID ANNOYING ERROR MESSAGE WHEN SETTING 'content':
	//// Error dojoUnique36_widget running custom onLoad code
	//fileManagerDialog._onLoadHandler = function(){
	//	this.isLoaded = true;
	//	try{
	//		this.onLoad.call(this);
	//	}catch(e){
	//		console.error('Error '+this.widgetId+' running custom onLoad code');
	//	}
	//};
*/

	// SET FILE MANAGER DIALOGUE CONTENT TO PANE NODE
	//fileManagerDialog.attr('content', paneNode);
	fileManagerDialog.content = paneNode;
	//fileManagerDialog.containerNode.appendChild(paneNode);


	//// ****************
	//// *** SOURCES ****
	//// ****************
	//var sources = this.getSources();
	//for ( var i = 0; i < sources.length; i++ )
	//{
	//	//if ( i > 0 )	{	break 	};
	//
	//	var source = sources[i];
	//
	//	// CREATE TITLE 
	//	var titlePaneNode = document.createElement('div');
	//	paneNode.appendChild(titlePaneNode);
	//
	//
	//	// CREATE TITLE PANE			
	//	var titlePane = new dijit.TitlePane(
	//		{
	//			title: "Source: " + sources[i].name,
	//			duration : 100, 
	//			open: false
	//		},
	//		titlePaneNode
	//	);
	//	
	//	var titleNodeId = dojo.dnd.getUniqueId();
	//	titlePane.titleNode.setAttribute( 'id', titleNodeId );
	//
	//	
	//	// ADD TOOLTIP
	//	
	//	//var tooltipLabel = "<span style='background-color: yellow'>" + sources[i].description + "</span>";
	//	var tooltipLabel = sources[i].description;
	//	var tooltip = new dijit.Tooltip(
	//		{
	//			label: tooltipLabel,
	//			connectId:[titleNodeId],
	//			position : ['above']
	//		}
	//	);
	//
	//	// SET URL
	//	var url = Agua.cgiUrl;
	//	url += "project.cgi?";
	//
	//	var query = "mode=fileSystem";
	//	query += "&username=" + Agua.cookie('username');
	//	query += "&sessionId=" + Agua.cookie('sessionId');
	//
	//
	//
	//	// SET LOCATION TO BE USED INSTEAD OF FILEROOT
	//	// ASSUMES ALL SOURCES ARE AT LEAST TWO FOLDERS DEEP IN THE FILESYSTEM
	//	
	//	/// NEXT IF NO MATCH
	//	if ( ! source.location.match(/^(.+)\/[^\/]+$/) )	continue;
	//		
	//	var sourceBase = source.location.match(/^(.+)\/[^\/]+$/)[1];
	//	var sourceFolder = source.location.match(/^.+\/([^\/]+)$/)[1];
	//	sourceBase = sourceBase.replace(/^\s+/,'');
	//	sourceBase = sourceBase.replace(/\s+$/,'');
	//
	//	query += "&location=" + sourceBase;
	//	
	//	// CREATE STORE
	//	var sourceStore = new dojox.data.FileStore(
	//		{
	//			url: url + query,
	//			pathAsQueryParam: true
	//		}
	//	);
	//
	//	// NOT NECCESSARY -- REMOVE LATER???
	//	////// SET FILE STORE path TO project
	//	////sourceStore.preamble = function()
	//	////{
	//	////	this.store.path = this.arguments[0].path;                        
	//	////};
	//
	//	// GENERATE NEW FileDrag OBJECT
	//	var sourceSelectorId = dojo.dnd.getUniqueId();
	//	var sourceSelector = new plugins.workflow.FileSelector(
	//		{
	//			id: sourceSelectorId,
	//			style: "height: auto; width: 100%; minHeight: 50px;",
	//			store: sourceStore,
	//			callback: this.callback,
	//			parentWidget: this,
	//			workflowObject: this.workflowObject
	//			
	//		}
	//	);
	//	
	//	// SET PATH FOR THIS PROJECT
	//	sourceSelector.path = sourceFolder;                    
	//	
	//	// SET LOCATION FOR USE IN CALLBACK WHEN FILE IS CLICKED
	//	// full path = location + file
	//	sourceSelector.location = sourceBase;
	//	//sourceSelector.attr('location', source.location);
	//	
	//	// START UP FileDrag
	//	sourceSelector.startup();
	//	
	//	// ADD sourceSelector TO TITLE PANE
	//	titlePane.containerNode.appendChild(sourceSelector.domNode);
	//
	//} // for loop on sources.length		
	//

	// POP UP DIALOGUE WINDOW
	fileManagerDialog.startup();
	//fileManagerDialog.show();

	// SHOW CONTENT
	// LATER: FIX THIS ERROR 'Error undefined running custom onLoad code'
	fileManagerDialog.set('content', paneNode);

	//fileManagerDialog.domNode.style.left='100px';
	fileManagerDialog.domNode.style.top='100px';


	// SET CSS CLASS
	dojo.addClass(fileManagerDialog.domNode, 'fileManager');

/*
	// SET DIALOGUE WRAPPER CLASS
	var wrapperId = dialogId + "_underlay_wrapper";
	var wrapperNode = dojo.byId(wrapperId);
	wrapperNode.setAttribute('class', 'fileManagerUnderlay dijitDialogUnderlayWrapper centered_underlay');
	//wrapperNode.setAttribute('class', 'fileManagerUnderlay');
*/
	// SET fileManagerDialog
	this.fileManagerDialog = fileManagerDialog;

}




}); 
