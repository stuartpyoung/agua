
dojo.provide("plugins.files.FileUpload");

// DISPLAY A FILE UPLOAD DIALOGUE AND ALLOW THE USER
// TO BROWSE FILES, SELECT AND THEN UPLOAD THEM

// DIALOGUE
dojo.require("dojox.widget.Dialog");

// FILE UPLOAD
dojo.require("dijit.form.Button");	
dojo.require("dojox.form.FileInput");	
dojo.require("dojox.form.FileInputAuto");	


dojo.require("dojox.form.FileUploader");
dojo.require("dijit.form.Button");
dojo.require("dijit.form.Textarea");
dojo.require("dojo.parser");

dojo.declare( "plugins.files.FileUpload",
	dojox.widget.RollingList,
{
	// FILE PATH TO UPLOAD FILE TO
	path : '',

	// STORE DND SOURCE ID AND DND TARGET ID
	sourceId: '',

	// DIALOGUE TO DISPLAY FILE MANAGER
	fileUploadDialog : null,

	// FILE UPLOADER OBJECT
	fileUploader : null,

	// ID FOR THIS FILE MANAGER DIALOG PANE
	dialogId : null,

	// UPLOAD BUTTON AND FILE TO UPLOAD TEXTBOX
	uploadButton : null,
	fileToUpload : null,

	// callback FUNCTION AND DATA FROM OBJECT THAT GENERATED THE FileUpload
	callback : null,

	preamble: function(){

	},

 	// CONSTRUCTOR	
	constructor : function(args) {

		this.callback = args.callback;

		// SET PANE ID
		this.path = args.path;

        // LOAD SORIA AND FILEPICKER CSS
        this.loadCSS();

		// LOAD APPLICATIONS DND SOURCE AND TARGET INTO FRAMEWORK
		this.showUploadDialog();
	},


	showUploadDialog: function ()
	{	

		var dialogId = dojo.dnd.getUniqueId();		
		this.dialogId = dialogId;

		var path = this.path;

		//// CREATE PANE NODE
		var containerNode = document.createElement('div');
		var containerNodeId = dojo.dnd.getUniqueId();
		containerNode.id = containerNodeId;

		var paneNode = document.createElement('div');
		var paneNodeId = dojo.dnd.getUniqueId();
		paneNode.id = paneNodeId;
		document.body.appendChild(paneNode);

		dojo.require("dojox.layout.ContentPane");
		var paneWidget = new dojox.layout.ContentPane(
			{
				executeScripts: true
			},
			paneNode
		);

		////using dojo.connect instead of addOnLoad worked
		////Submitted by dwrobertson on Tue, 03/04/2008 - 20:03.
		////Thank you for the fix. Using dojo.connect(dijit.byId('testPane'),'onLoad',init);
		////instead of dojo.addOnLoad(init) worked for the dojox.layout.ContentPane.
		////

//		var handle = paneWidget.setHref('plugins/project/FileUpload/fileUpload.html');
//		handle.addOnLoad(function(){
//


		// REPLACE HREF WITH PROGRAMMATIC CONSTRUCTION OF FILE UPLOAD
		//var handle = paneWidget.setHref('plugins/project/FileUpload/fileUpload.html');

		var uploadObject = this;
		if ( uploadObject.table == null )
		{
			var table = document.createElement('table');
			this.table = table;
			//document.body.appendChild(table);
			//paneWidget.attr('content', table);
			paneWidget.domNode.appendChild(table);

			var row = document.createElement('row');
			table.appendChild(row);

			var td = document.createElement('td');
			row.appendChild(td);
			td.setAttribute('class', 'leftCol');

			// ADD SELECT FILES BUTTON
			var selectButtonId = dojo.dnd.getUniqueId();
			var selectButton = new dijit.form.Button(
				{
					label : "Select Files..."
					,
					id: selectButtonId
				}
			);
			selectButton.attr('class', "browse");
			td.appendChild(selectButton.domNode);
			selectButton.domNode.onmouseover = function(){};
			selectButton.domNode.onmouseout = function(){};


			dojo.connect(selectButton, "onClick", function(e)
			{
			});


			// SET SELECT BUTTON DOMNODE ID
			selectButton.domNode.setAttribute('id', selectButtonId);

			// SET this.selectButton
			uploadObject.selectButton = selectButton;

			// ADD UPLOAD BUTTON
			var uploadButton = new dijit.form.Button(
				{
					label : "upload Files..."
				}
			);
			uploadButton.attr('class', "uploadBtn");
			td.appendChild(uploadButton.domNode);

			// SET this.uploadButton
			uploadObject.uploadButton = uploadButton;

			// SET UPLOAD BUTTON ONCLICK
			var fileUploader = this;
			dojo.connect(uploadButton, "onClick", function(e)
			{
				doUpload();
			});

			// SECOND ROW WITH FILE INPUT
			var row2 = document.createElement('row');
			table.appendChild(row2);

			var td2 = document.createElement('td');
			row2.appendChild(td2);
			td2.setAttribute('class', 'leftCol');

			// TEXT AREA TO DISPLAY UPLOAD FILE NAME
			var textAreaContainer = document.createElement('div');
			td2.appendChild(textAreaContainer);
			var fileToUpload = new dijit.form.Textarea(
				{
					value: "File to upload goes here"
				},
				textAreaContainer
			);

			// SET this.fileToUpload
			uploadObject.fileToUpload = fileToUpload;

			var td3 = document.createElement('td');
			row2.appendChild(td3);
			td3.setAttribute('class', 'rgtCol');


			//var uploadUrl = "dojo-1.5.0/dojox/form/resources/UploadFile.php";
			var uploadUrl = Agua.cgiUrl + "upload.cgi";

			//dojo.byId("uploadedFiles").value = "";
			this.fileToUpload.value = "";
		}


		var f0 = new dojox.form.FileUploader({
			button: this.selectButton, 
			degradable:true,
			uploadUrl:uploadUrl, 
			uploadOnChange:false, 
			selectMultipleFiles: false,
			//fileMask:fileMask,
			isDebug:true
		});

		this.fileUploader = f0;



		// OVERRIDE upload FUNCTION TO ADD path TO POST VARIABLES
		f0.upload = function(){


			var childNodes = this._formNode.childNodes;
			var notFound = 1;
			if ( childNodes )
			{
				for ( var i = 0; i < childNodes.length; i++ )
				{
					if ( childNodes[i].id == "hiddenUsername" )
					{
						notFound = 0;
						break;
					}
				}
			}

			// APPEND HIDDEN NODES TO FORM
			if ( notFound )
			{
				var userNode = document.createElement('input');
				userNode.type = "hidden";
				userNode.name = "username";
				userNode.value = Agua.cookie('username');
				userNode.id = "hiddenUsername";
				this._formNode.appendChild(userNode);

				var sessionNode = document.createElement('input');
				sessionNode.type = "hidden";
				sessionNode.name = "sessionId";
				sessionNode.value = Agua.cookie('sessionId');
				sessionNode.id = "hiddenSessionId";
				this._formNode.appendChild(sessionNode);

				var pathNode = document.createElement('input');
				pathNode.type = "hidden";
				pathNode.name = "path";
				pathNode.value = path;
				pathNode.id = "hiddenPath";
				this._formNode.appendChild(pathNode);
			}


			// summary
			//	Tell form to upload
			dojo.io.iframe.send({
				//url: this.uploadUrl + "?username=admin&sessionId=1228791394.7868.158",
				url: this.uploadUrl,
				form: this._formNode,
				timeout: 20000,
				handleAs: "json",
				handle: dojo.hitch(this,function(data,ioArgs,widgetRef){
					this.onComplete(this.selectMultipleFiles?data:[data]);								 
				})	
			});

		},



		doUpload = function(){


			uploadObject.fileToUpload.innerHTML = "uploading...";
			f0.upload();
		}

		dojo.connect(f0, "onChange", function(data){
			dojo.forEach(data,
			function(d)
			{
				dojo.byId("fileToUpload").value = d.name+" "+Math.ceil(d.size*.001)+"kb \n";
			});
		});

		dojo.connect(f0, "onProgress", function(data){
			dojo.byId("fileToUpload").value = "";
			dojo.forEach(data, function(d){
				this.fileToUpload.value += "("+d.percent+"%) "+d.name+" \n";

			});
		});

		dojo.connect(f0, "onComplete", function(data){
			dojo.forEach(data, function(d){
				//dojo.byId("uploadedFiles").value += d.file+" \n";
				dojo.byId("rgtCol").innerHTML += '<img src="'+d.file+'" />';
			});
		});

		Destroy = function(){
			f0.destroyAll();
		}



		//handle.addOnUnload(function(){
		//	cleanUp();
		//});




		// ATTACH HTML TO CONTAINER NODE OF FILE UPLOAD DIALOG

		//dojo.require("dojo.parser");
		//dojo.parser.parse(paneNode.domNode);

		// CREATE DIALOGUE
		var fileUploadDialog = new dojox.widget.Dialog ({
			id: dialogId,
			dimensions: [800,"auto"],
			draggable: true,
			align: 'left',
			sizeDuration: 200,
			sizeMethod:       "combine",
			viewportPadding : "125",
			showTitle: "true",
			title: "File Upload",
			executeScripts: true

		}, "fileUploadDialog" );		
		fileUploadDialog.attr('widgetId', dialogId + "_widget");


		// DEBUG

		// SET CLASS FOR STYLE OF INSERTED PANE
		dojo.attr(dojo.byId(dialogId), 'class', 'fileUpload dijitDialog');
		dojo.attr(fileUploadDialog.titleBar, 'class', 'fileUploadTitle');
		//dojo.marginBox(fileUploadDialog.titleBar).h = 30px;


		// SHOW CONTENT
		// LATER: FIX T	HIS ERROR 'Error undefined running custom onLoad code'
		//fileUploadDialog.setContent(paneNode);

		// POP UP DIALOGUE WINDOW
		fileUploadDialog.startup();
		fileUploadDialog.show();

		//fileUploadDialog.domNode.style.left='100px';
		fileUploadDialog.domNode.style.top='100px';


		//document.body.appendChild(paneNode);
		//fileUploadDialog.setContent(paneNode);


		fileUploadDialog.dojoxDialogWrapper.appendChild(containerNode);
		//document.body.appendChild(containerNode);

		containerNode.appendChild(paneNode);


			//var f0 = new dojox.form.FileUploader({
			//	button:dijit.byId("btn0"), 
			//	degradable:true,
			//	uploadUrl:uploadUrl, 
			//	uploadOnChange:false, 
			//	selectMultipleFiles: false,
			//	fileMask:fileMask,
			//	isDebug:true
			//});
			//


		//// SET DIALOGUE WRAPPER CLASS
		//var wrapperId = dialogId + "_underlay_wrapper";
		//var wrapperNode = dojo.byId(wrapperId);
		//wrapperNode.setAttribute('class', 'fileUploadUnderlay dijitDialogUnderlayWrapper centered_underlay');
		//wrapperNode.setAttribute('class', 'fileUploadUnderlay');


		//dojo.require("dojo.parser");
		//dojo.parser.parse(paneNode);


		//fileUploadDialog.onLoad(function(){

		//	dojo.require("dojox.form.FileUploader");
		//	dojo.require("dijit.form.Button"); 
		//	dojo.require("dojo.parser");
		//
		//	var uploadUrl = "../dojo-1.5.0/dojox/form/resources/UploadFile.php";
		//
		//	var rmFiles = "";
		//	var fileMask = [
		//		["Jpeg File", 	"*.jpg;*.jpeg"],
		//		["GIF File", 	"*.gif"],
		//		["PNG File", 	"*.png"],
		//		["All Images", 	"*.jpg;*.jpeg;*.gif;*.png"]
		//	];
		//	// For testing 1D array masks:
		//	// var fileMask = 	["All Images", 	"*.jpg;*.jpeg;*.gif;*.png"];
		//	// var fileMask = 	["PNG File", 	"*.png"];
		//
		//
		//	dojo.byId("uploadedFiles").value = "";
		//	dojo.byId("fileToUpload").value = "";
		//
		//
		//	var f0 = new dojox.form.FileUploader({
		//		button:dijit.byId("btn0"), 
		//		degradable:true,
		//		uploadUrl:uploadUrl, 
		//		uploadOnChange:false, 
		//		selectMultipleFiles: false,
		//		fileMask:fileMask,
		//		isDebug:true
		//	});
		//
		//	doUpload = function(){
		//		dojo.byId("fileToUpload").innerHTML = "uploading...";
		//		f0.upload();
		//	}
		//
		//	dojo.connect(f0, "onChange", function(data){
		//		dojo.forEach(data, function(d){
		//
		//				dojo.byId("fileToUpload").value = d.name+" "+Math.ceil(d.size*.001)+"kb \n";
		//		});
		//	});
		//
		//	dojo.connect(f0, "onProgress", function(data){
		//		dojo.byId("fileToUpload").value = "";
		//		dojo.forEach(data, function(d){
		//			dojo.byId("fileToUpload").value += "("+d.percent+"%) "+d.name+" \n";
		//			
		//		});
		//	});
		//
		//	dojo.connect(f0, "onComplete", function(data){
		//		dojo.forEach(data, function(d){
		//			dojo.byId("uploadedFiles").value += d.file+" \n";
		//			dojo.byId("rgtCol").innerHTML += '<img src="'+d.file+'" />';
		//			rmFiles+=d.file+";";
		//		});
		//	});
		//
		//	Destroy = function(){
		//		f0.destroyAll();
		//	}
		//
		//	cleanUp = function()
		//	{
		//		dojo.byId("rgtCol").innerHTML = "";
		//		dojo.byId("uploadedFiles").value = "";
		//		dojo.xhrGet({
		//			url:uploadUrl,
		//			handleAs:"text",
		//			content:{
		//				rmFiles:rmFiles
		//			}
		//		});
		//		rmFiles = "";
		//	}
		//
		//
		//	dojo.addOnUnload(function(){
		//		cleanUp();
		//	});
		//	
		////});
		//

		// SET fileUploadDialog
		this.fileUploadDialog = fileUploadDialog;
	},


    loadCSS : function ()
    {        

        // THIS CSS FILE PROVIDES THE ICONS AND FORMATTING FOR INDIVIDUAL FILES/DIRECTORIES
		var cssFile3 = "plugins/project/FileUpload/fileUpload.css";
		var cssNode = document.createElement('link');
		cssNode.type = 'text/css';
		cssNode.rel = 'stylesheet';
		cssNode.href = cssFile3;
		cssNode.media = 'screen';
		cssNode.title = 'loginCSS';
		cssNode.id = "widgetStyle";
		document.getElementsByTagName("head")[0].appendChild(cssNode);
    },


	getSyncJson : function (url)
	{

		var json;
		dojo.xhrGet(
			{
				url: url,
				//handleAs: "json",
				handleAs: "json-comment-optional",
				sync: true,
				handle: function(response){
					json = response;
				}
			}
		);

		return json;
	}


}); 
