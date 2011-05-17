
dojo.provide("plugins.files.FileMenu");

// WIDGET PARSER
dojo.require("dojo.parser");

// INHERITS
dojo.require("plugins.core.Common");

// HAS A
dojo.require("plugins.menu.Menu");
dojo.require("plugins.core.InputDialog");
dojo.require("plugins.core.InteractiveDialog");
dojo.require("plugins.core.ConfirmDialog");
dojo.require("dojox.form.FileInputAuto");


dojo.declare(
    "plugins.files.FileMenu",
	[ dijit._Widget, dijit._Templated, plugins.core.Common ],
{
	/////}}

//Path to the template of this widget. 
templatePath: dojo.moduleUrl("plugins", "files/templates/filemenu.html"),

// Calls dijit._Templated.widgetsInTemplate
widgetsInTemplate : true,

//addingApp STATE
addingApp : false,

// OR USE @import IN HTML TEMPLATE
cssFiles : [ dojo.moduleUrl("plugins.files") + "/css/filemenu.css" ],

// PARENT WIDGET
parentWidget : null,

constructor : function(args) {
	// GET INFO FROM ARGS
	//this.parentWidget = args.parentWidget;

	// LOAD SORIA AND FILEPICKER CSS
	this.loadCSS();		
},

postCreate : function() {

	// SET INPUT DIALOG
	this.setInputDialog();

	// SET INTERACTIVE DIALOG
	this.setInteractiveDialog();

	// SET CONFIRM DIALOG
	this.setConfirmDialog();

	// SET LABEL
	this.setTitle("File Menu");

	// CONNECT SHORTKEYS FOR MENU
	this.setMenu();

	// SET THE UPLOAD OBJECT
	this.setUploader();

	// DO INHERITED STARTUP
	this.startup();
},


startup : function () {

	// COMPLETE CONSTRUCTION OF OBJECT
	this.inherited(arguments);	 

},

// CREATE BLIND FILE UPLOAD WIDGET IN uploadNode MENU ITEM
// NB: MUST USE CSS TO HIDE AND STUNT UNNECESSARY NODES: E.G., fakeNodeHolder
setUploader : function () {

	// CREATE BLIND FILE INPUT INSIDE uploadNode MENU ITEM
	var inputNode = document.createElement('input');
	inputNode.setAttribute('size', 5);
	this.uploadNode.containerNode.appendChild(inputNode);

	var uploadUrl = Agua.cgiUrl + "upload.cgi";
	this.blindUpload = new dojox.form.FileInputBlind(
		{
			url : uploadUrl
		},
		inputNode
	);

	// STOP MENU CLOSE WHEN ONCLICK EVENT REGISTERS FIRST AT MENU ITEM
	this.uploadNode._onClick = function(event)
	{

		//dojo.stopEvent(event);
		return;
	}

	//// **** DEBUG ****
	//
	//// STOP CLOSE
	//this.uploadNode._blur = function(event)
	//{
	//
	//	dojo.stopEvent(event);
	//	return;
	//}
	//

	//// STOP CLOSE
	//this.menu._onBlur = function(event)
	//{
	//
	//	dojo.stopEvent(event);
	//	return;
	//}


	dojo.connect(this.uploadNode, "onClick", dojo.hitch(this, function(event)
		{

			event.stopPropagation();
		}
	));


	// TRIGGER THE CHAIN OF EVENTS TO UPLOAD A FILE IN THE BACKGROUND
	var thisObject = this;
	this.blindUpload._sendFile = function(event)
	{

		// TOGGLE SENDING STATE
		if(this._sent || this._sending || !this.fileInput.value){ return; }
		this._sending = true;

		dojo.fadeIn({ node: this.overlay, duration:this.duration }).play();

		var _newForm; 
		if(dojo.isIE){
			// just to reiterate, IE is a steaming pile of code. 
			_newForm = document.createElement('<form enctype="multipart/form-data" method="post">');
			_newForm.encoding = "multipart/form-data";

		}else{
			// this is how all other sane browsers do it
			_newForm = document.createElement('form');
			_newForm.setAttribute("enctype","multipart/form-data");
		}
		_newForm.appendChild(this.fileInput);
		dojo.body().appendChild(_newForm);



		// GET THE PATH
		var path = thisObject.getPath();

		// SET HIDDEN VALUES TO GENERATE HIDDEN NODES IN FORM
		var hiddenValues = {
			username: Agua.cookie('username'),
			sessionId : Agua.cookie('sessionId'),
			path: path
		};		

		// CREATE HIDDEN NODES
		for ( var key in hiddenValues )
		{
			var node = document.createElement('input');
			node.type = "hidden";
			node.name = key;
			node.value = hiddenValues[key];
			node.id = "hidden" + key;
			_newForm.appendChild(node);
		}


		dojo.io.iframe.send({
			url: this.url,
			form: _newForm,
			timeout : 1000,
			handleAs: "json",
			//handle: dojo.hitch(thisObject,"onUploadComplete")
			handle: "_handleSend"
		});

	};


	// summary: The callback to toggle the progressbar, and
	// fire the user-defined callback
	this.blindUpload._handleSend = function(data,ioArgs)
	{	

		// DO onUploadComplete
		thisObject.onUploadComplete();

		// innerHTML throws errors in IE! so use DOM manipulation instead
		this.overlay.removeChild(this.overlay.firstChild);

		this._sent = true;
		this._sending = false;
		dojo.style(this.overlay,{
			opacity:0,
			border:"none",
			background:"none"
		}); 

		this.overlay.style.backgroundImage = "none";
		this.fileInput.style.display = "none";
		this.fakeNodeHolder.style.display = "none";
		dojo.fadeIn({ node:this.overlay, duration:this.duration }).play(250);

		this.disconnect(this._blurListener);
		this.disconnect(this._focusListener);

		//remove the form used to send the request
		dojo.body().removeChild(ioArgs.args.form);
		//this.fileInput = null;

		this.onComplete(data,ioArgs,this);
	};


	// fired when an upload has finished. 
	this.blindUpload.onComplete = function(data,ioArgs,widgetRef)
	{
		// data: the raw data found in the first [TEXTAREA] tag of the post url
		// ioArgs: the dojo.Deferred data being passed from the handle: callback
		// widgetRef: this widget pointer, so you can set this.overlay to a completed/error message easily


	}
},	//	setUploader


// RETURN THE FILE PATH OF THE FOCUSED GROUP DRAG PANE
getPath : function () {
	var groupDragPane = dijit.getEnclosingWidget(this.menu.currentTarget);

	if ( groupDragPane == null || ! groupDragPane )	return;

	return groupDragPane.path;
},


// SET THE MENU TITLE
setTitle : function (title) {
// NO TITLE - DO NOTHING


	//this.titleNode.containerNode.innerHTML = title;
},

// CONNECT SHORTKEYS FOR MENU
setMenu : function () {

	// NOTE: USE accelKey IN DOJO 1.3 ONWARDS
	dojo.connect(this.menu, "onKeyPress", dojo.hitch(this, function(event)
	{
		var key = event.charOrCode;
		if ( this.altOn == true )
		{
			switch (key)
			{
				case "n" : this.newFolder(); break;
				case "m" : this.rename(); break;
				case "l" : this.deleteNode(); break;
				case "o" : this.openWorkflow(); break;
				case "u" : this.upload(); break;
				case "w" : this.download(); break;
			}
		}
		event.stopPropagation();
	}));

	// SET ALT KEY ON/OFF
	dojo.connect(this.menu, "onKeyDown", dojo.hitch(this, function(event){
		var keycode = event.keyCode;
		if ( keycode == 18 )	this.altOn = true;
	}));
	dojo.connect(this.menu, "onKeyUp", dojo.hitch(this, function(event){
		var keycode = event.keyCode;
		if ( keycode == 18 )	this.altOn = false;
	}));

	this.menu.onCancel = function(event)
	{
	}

},

// RELOAD THE GROUP DRAG PANE ON UPLOAD COMPLETE
onUploadComplete : function () {
	var groupDragPane = dijit.getEnclosingWidget(this.menu.currentTarget);
return;

	// RELOAD THE PROJECTS TAB
	setTimeout(function(thisObj) { groupDragPane.reload(); }, 1000, this);
},


// BIND THE MENU TO A NODE
bind : function (node) {

	if ( node == null )
	{

	}
	return this.menu.bindDomNode(node);	
},


// openWorkflow THE INPUTS AND OUTPUTS OF THIS STAGE TO THE PARAMETER VALUES
// OF THE PRECEDING STAGE
openWorkflow : function () {

},

// RENAME FILE	
rename : function () {
	var oldFilename = this.menu.currentTarget.innerHTML;

	var groupDragPane = dijit.getEnclosingWidget(this.menu.currentTarget);

	// SET TITLE AND MESSAGE
	var title = "Rename file '" + oldFilename + "'";
	var message = "Please enter new name";

	// CALLBACKS
	var cancelCallback = function (){
	};
	var enterCallback = dojo.hitch(this, function (newFilename)
		{

			// SANITY CHECK
			if ( newFilename == null )	return;
			newFilename = newFilename.replace(/\s+/, '');
			if ( newFilename == '' )	return;

			// CHECK IF NAME EXISTS ALREADY
			if ( this.inFiles(groupDragPane, newFilename) == true )
			{
				return;
			}

			this.menu.currentTarget.innerHTML = newFilename;
		}
	);		

	// SHOW THE DIALOG
	this.loadInputDialog(title, message, enterCallback, cancelCallback);
},

// RETURN TRUE IF FILE NAME EXISTS IN THIS GROUP DRAG PANE, ELSE RETURN FALSE
inFiles : function (groupDragPane, newFilename) {
	if ( newFilename == null )	return;

	var childNodes = groupDragPane._dragSource.getAllNodes();

	for ( var i = 0; i < childNodes.length; i++ )
	{
		if ( newFilename == childNodes[i].innerHTML )	return true;
	}

	return false;
},


// CREATE A NEW FOLDER
newFolder : function () {
	if ( this.menu.currentTarget == null )	return;

	var groupDragPane = dijit.getEnclosingWidget(this.menu.currentTarget);

	// SET TITLE AND MESSAGE
	var title = "New Folder";
	var message = "Please enter folder name";

	// CALLBACKS
	var cancelCallback = function (){
	};
	var enterCallback = dojo.hitch(this, function (newFoldername)
		{

			// SANITY CHECK
			if ( newFoldername == null )	return;
			newFoldername = newFoldername.replace(/\s+/, '');
			if ( newFoldername == '' )	return;

			// CHECK IF NAME EXISTS ALREADY
			if ( this.inFiles(groupDragPane, newFoldername) == true )
			{
				return;
			}

			// INSERT NEW CHILD
			var item = new Object;
			item.name = newFoldername;
			item.type = ["folder"];
			item.directory = true;
			item.parentPath = groupDragPane.path;
			item.path = newFoldername;
			groupDragPane.items.push(item);

			// GENERATE CHILD
			var newChild = groupDragPane.parentWidget._getMenuItemForItem(item, groupDragPane);

			// INSERT CHILD
			groupDragPane._dragSource.insertNodes(false, [ newChild ]);


			// ADD 'directory' CLASS TO NEW CHILD
			var allNodes = groupDragPane._dragSource.getAllNodes();
			var node = allNodes[allNodes.length - 1];
			node.setAttribute("dndType", "file");
			dojo.addClass(node, "directory");

			// ADD item ATTRIBUTE ON NODE
			node.item = item;
			node.item._S = newChild.store;

			// ADD MENU
			var dynamicMenu = groupDragPane.createMenu("folder");

			// BIND THE MENU TO THE DND NODE
			dynamicMenu.bindDomNode(node);

			// CONNECT ONCLICK
			dojo.connect(node, "onclick", groupDragPane, "onclickHandler");

			// CREATE FOLDER ON SERVER	
			var url = Agua.cgiUrl + "project.cgi?";
			var folderPath = groupDragPane.path + "/" + newFoldername;
			var query = "mode=newFolder";
			query += "&sessionId=" + Agua.cookie('sessionId');
			query += "&username=" + Agua.cookie('username');
			query += "&folderpath=" + folderPath;

			dojo.xhrGet(
				{
					url: url + query,
					handleAs: "json",
						//handleAs: "json-comment-optional",
					sync: false,
					handle: function(response){
						if ( response.error )
							{
						}
						else
						{
						}
					}
				}
			);
		}
	);		

	// SHOW THE DIALOG
	this.loadInputDialog(title, message, enterCallback, cancelCallback);

},	//	newFolder



// deleteNode SELECTED NODE
deleteNode : function () {
	if ( this.menu.currentTarget == null )	return;

	// GET THE PROJECT WIDGET
	var filename = this.menu.currentTarget.item.name;

	var isDirectory = this.menu.currentTarget.item.directory;
	var type = "file";
	if ( isDirectory == true )	type = "folder";

	var groupDragPane = dijit.getEnclosingWidget(this.menu.currentTarget);

	// CALLBACKS
	var noCallback = function (){
	};
	var yesCallback = dojo.hitch(this, function ()
		{
			groupDragPane._dragSource.deleteSelectedNodes();

			var item = this.menu.currentTarget.item;
			var file = item.parentPath + "/" + item.path;

			var url = Agua.cgiUrl + "project.cgi?";
			var query = "mode=removeFile";
			query += "&sessionId=" + Agua.cookie('sessionId');
			query += "&username=" + Agua.cookie('username');
			query += "&file=" + file;
			//var query = "mode=removeFile&sessionId=1228791394.7868.158&username=admin&file=" + file;

			dojo.xhrGet(
				{
					url: url + query,
					handleAs: "json",
						//handleAs: "json-comment-optional",
					sync: false,
					handle: function(response){
						if ( response.error )
						{
							// status: initiated, ongoing, completed
							// error: file not found, removing, file being copied
							fileCopyError = response.error;
						}
						if ( ! response.error )
						{
							var parentNode = this.menu.currentTarget.parentNode;
							parentNode.removeChild(this.menu.currentTarget);
						}
					}
				}
			);
		}
	);
			// SET TITLE AND MESSAGE
	var title = "Delete " + type + ": '" + filename + "'?";
	var message = "All its data will be destroyed";

	// SHOW THE DIALOG
	this.loadConfirmDialog(title, message, yesCallback, noCallback);

},	//	deleteNode

// DOWNLOAD FILE FROM FOLDER
download : function () {

},



openFileDownload : function (filepath) {

	var query = "?username=" + Agua.cookie('username');
	query += "&sessionId=" + Agua.cookie('sessionId');
	query += "&filepath=" + filepath;

	var url = Agua.cgiUrl + "download.cgi";

	var args = {
		method: "GET",
		url: url + query,
		//content: {},
		handleAs: "html",
		timeout: 10000
		//load: dojo.hitch(this, "onDownloadComplete"),
		//error: dojo.hitch(this, "onDownloadError")
	};

	// do an IFrame request to download the csv file.
	dojo.io.iframe.send(args);
},



commitRename : function (item, newName) {

	var url = Agua.cgiUrl + "project.cgi?";
	var file = item.parentPath + "/" + item.path;

	var query = "mode=commitRename";
	query += "&sessionId=" + Agua.cookie('sessionId');
	query += "&username=" + Agua.cookie('username');
	query += "&file=" + file;
	query += "&newName=" + newName;

	dojo.xhrGet(
		{
			url: url + query,
			handleAs: "json",
				//handleAs: "json-comment-optional",
			sync: false,
			handle: function(response){
				if ( response.error )
				{
					// status: initiated, ongoing, completed
					// error: file not found, removing, file being copied
					fileCopyError = response.error;
				}
				if ( ! response.error )
				{
					var parentNode = this.menu.currentTarget.parentNode;
					parentNode.removeChild(this.menu.currentTarget);
				}
			}
		}
	);
},

setConfirmDialog : function () {
	var yesCallback = function (){};
	var noCallback = function (){};
	var title = "Dialog title";
	var message = "Dialog message";

	this.confirmDialog = new plugins.core.ConfirmDialog(
		{
			title 				:	title,
			message 			:	message,
			parentWidget 		:	this,
			yesCallback 		:	yesCallback,
			noCallback 			:	noCallback
		}			
	);
},

loadConfirmDialog : function (title, message, yesCallback, noCallback) {

	this.confirmDialog.load(
		{
			title 				:	title,
			message 			:	message,
			yesCallback 		:	yesCallback,
			noCallback 			:	noCallback
		}			
	);
},


setInputDialog : function () {
	var enterCallback = function (){};
	var cancelCallback = function (){};
	var title = "";
	var message = "";

	this.inputDialog = new plugins.core.InputDialog(
		{
			title 				:	title,
			message 			:	message,
			parentWidget 		:	this,
			enterCallback 		:	enterCallback,
			cancelCallback 		:	cancelCallback
		}			
	);
},

loadInputDialog : function (title, message, enterCallback, cancelCallback) {

	this.inputDialog.load(
		{
			title 				:	title,
			message 			:	message,
			enterCallback 		:	enterCallback,
			cancelCallback 		:	cancelCallback
		}			
	);
},

setInteractiveDialog : function () {
	var enterCallback = function (){};
	var cancelCallback = function (){};
	var title = "";
	var message = "";

	this.interactiveDialog = new plugins.core.InteractiveDialog(
		{
			title 				:	title,
			message 			:	message,
			parentWidget 		:	this,
			enterCallback 		:	enterCallback,
			cancelCallback 		:	cancelCallback
		}			
	);
},

loadInteractiveDialog : function (title, message, enterCallback, cancelCallback, checkboxMessage) {

	this.interactiveDialog.load(
		{
			title 				:	title,
			message 			:	message,
			checkboxMessage 	:	checkboxMessage,
			enterCallback 		:	enterCallback,
			cancelCallback 		:	cancelCallback
		}			
	);
},

hide : function () {

	dojo.style(this.containerNode, {
		opacity: 0,
		overflow: "hidden"
	});
},

show : function () {

	dojo.style(this.containerNode, {
		opacity: 1,
		overflow: "visible"
	});
},

disable : function () {

	this.menu.enabled = false;
},

enable : function () {

	this.menu.enabled = true;
} 





}); // plugins.files.FileMenu
