
dojo.provide("plugins.upload.FileUp");

// DISPLAY A FILE UPLOAD DIALOGUE AND ALLOW THE USER
// TO BROWSE FILES, SELECT AND THEN UPLOAD THEM

// DIALOGUE
dojo.require("dojox.widget.Dialog");

// FILE UPLOAD
dojo.require("dijit.form.Button");	
dojo.require("dojox.form.FileInput");	
dojo.require("dojox.form.FileInputAuto");	

dojo.declare( "plugins.upload.FileUp",
	dojox.widget.RollingList,
{
	// FILE PATH TO UPLOAD FILE TO
	path : '',

	// STORE DND SOURCE ID AND DND TARGET ID
	sourceId: '',

	// DIALOGUE TO DISPLAY FILE MANAGER
	fileUploadDialog : null,

	// ID FOR THIS FILE MANAGER DIALOG PANE
	dialogId : null,

	// THE UPLOAD dojox.layout.ContentPane WIDGET
	paneWidget : null,

	// FILE UPLOADER
	fileUploader : null,

	// FILE TO BE UPLOADED
	uploadedFile : null,

	// callback FUNCTION AND DATA FROM OBJECT THAT GENERATED THE FileUpload
	onCompleteCallback : null,

	// PARENT WIDGET THAT HAS_A FileUp OBJECT
	parentWidget : null,


	preamble: function(){

	},


 	// CONSTRUCTOR	
	constructor : function(args) {

		// SET PARENT WIDGET
		this.parentWidget = args.parentWidget;

		// SET PANE ID
		this.path = args.path;

		// OPTION TABLE DATA NODE THE FILE IS BEING UPLOADED FOR
		this.node = args.node;

        // LOAD SORIA AND FILEPICKER CSS
        this.loadCSS();

		// LOAD APPLICATIONS DND SOURCE AND TARGET INTO FRAMEWORK
		this.showUploadDialog();
	},





	getFileInfo : function(filepath)
	{

		// SHOW LOADING...
		dojo.byId("fileUploadMessage").innerHTML = "Checking if file already exists ...";
		dojo.byId("fileUploadStatusIcon").setAttribute('class', 'fileUploadStatusIcon');

		// DO xhrPut TO OBTAIN FILE INFORMATION IF EXISTS ON FILE SYSTEM
		var fileInfo;

		// GET URL
		var url = Agua.cgiUrl + "workflow.cgi";

		// GENERATE QUERY JSON FOR THIS WORKFLOW IN THIS PROJECT
		var query = new Object;
		query.username = Agua.cookie('username');
		query.sessionId = Agua.cookie('sessionId');
		query.path = filepath;
		query.mode = "fileInfo";


		// SEND TO SERVER
		dojo.xhrPut(
			{
				url: url,
				contentType: "text",
				sync : true,
				handleAs: "json",
				//handleAs: "json-comment-filtered",
				putData: dojo.toJson(query),
				timeout: 30000,
				load: function(response, ioArgs) {
					fileInfo = response;
				},
				error: function(response, ioArgs) {
					return response;
				}
			}
		);	

		//fileInfo = {
		//	'exists'	: 	true,
		//	size		:	1000000,
		//	modified	:	'Tues 19 Aug 2009 21:00:01',
		//	directory	:	false
		//};

		// SHOW END LOADING
		dojo.byId("fileUploadMessage").innerHTML = "";
		dojo.byId("fileUploadStatusIcon").setAttribute('class', '');

		return fileInfo;
	},


	//destroyAll : function ()
	//{
	//	// DESTROY THE dojox.form.fileUploader OBJECT
	//	fileUploadObject.fileUploader.destroyAll();
	//	
	//	// HIDE THE FILE UPLOAD DIALOGUE
	//	this.fileUploadDialog.hide();
	//	
	//},


	confirmOverwrite : function (fileInfo, callback, data)
	{


		var fileUploadObject = this;


		// GENERATE THE CONFIRMATION DIALOGUE
		function raiseQueryDialog(title, question, callbackFn) {


			var ids = ["queryDialog", "yesButton", "noButton"];
			for ( var i = 0; i < ids.length; i++ )
			{
				if ( dijit.byId(ids[i]) )
				{
					dijit.byId(ids[i]).destroy();
				}
			}


			var node = document.createElement('div');
			document.body.appendChild(node);

			var confirmDialog = new dijit.Dialog({ id: 'queryDialog', title: title }, node);
			confirmDialog.attr('class', 'fileUploadConfirmationDialogue');
			dojo.addClass(confirmDialog.titleBar, 'fileUploadConfirmationDialogueTitle');

			// When either button is pressed, kill the dialog and call the callbackFn.
			var commonCallback = function(mouseEvent)
			{
				confirmDialog.hide();
				confirmDialog.destroyRecursive();
				if (mouseEvent.explicitOriginalTarget.id == 'yesButton')
				{
					callback(data);	
				}
				else
				{
					dojo.byId("fileToUpload").value = "";

					//// LATER: GET THROUGH fileUploadObject.fileUploader (i.e., f0)

					// DESTROY THE f0 OBJECT
					//callbacks.no();
					//fileUploadObject.fileUploader.destroyAll();
					//fileUploadObject.destroy();

					// HIDE THE FILE UPLOAD DIALOGUE
					fileUploadObject.fileUploadDialog.hide();

					return false;
				}
			};

			// CREATE CONTENT AND BUTTONS
			var questionDiv = document.createElement('div');
			questionDiv.innerHTML = question;
			var yesButton = new dijit.form.Button(
				{ label: 'Yes', id: 'yesButton', onClick: commonCallback }
			);
			yesButton.attr('class', 'fileUploadConfirmationButton');
			dojo.addClass(yesButton, 'fileUploadYesButton');

			var noButton = new dijit.form.Button(
				{ label: 'No', id: 'noButton', onClick: commonCallback }
			);
			noButton.attr('class', 'fileUploadConfirmationButton');
			dojo.addClass(noButton, 'fileUploadNoButton');

			// ADD CONTENT AND BUTTONS
			dojo.addClass(yesButton, 'fileUploadYesButton');
			yesButton.attr('class', 'fileUploadConfirmationButton');

			confirmDialog.containerNode.appendChild(questionDiv);
			confirmDialog.containerNode.appendChild(yesButton.domNode);
			confirmDialog.containerNode.appendChild(noButton.domNode);

			confirmDialog.show();
		}

		var type = "File";
		if ( fileInfo.directory )
		{
			type = "Directory";
		}
		var message = "<table class='fileUploadConfirmTable'><tr><td colspan='2' class='fileUploadConfirmTableIntro'>" + type + " already exists:</td></tr>";
		message += "<tr><td class='fileUploadConfirmTableHeader' colspan='2'> " + fileInfo.filepath + "</td></tr>";
		message += "<tr><td class='fileUploadConfirmTableItem'>Modified </td><td class='fileUploadConfirmTableValue'>" + fileInfo.modified + "</td></tr>";
		message += "<tr><td class='fileUploadConfirmTableItem'>Size     </td><td class='fileUploadConfirmTableValue'>" + String(fileInfo.size/1000) + " kb </td></tr>";
		message += "<tr><td colspan='2' class='fileUploadConfirmTableIntro'> Are you sure you want to overwrite this file? </td></tr>";
		message += "</table>";

		raiseQueryDialog("Confirm File Overwrite", message);
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

		dojo.require("dojox.layout.ContentPane");
		var paneWidget = new dojox.layout.ContentPane(
			{
				executeScripts: true
			},
			paneNode
		);

		// SET this.paneWidget
		this.paneWidget = paneWidget;


		// REMOVE TIDYING LATER
		//uploadBtn
		//btn0  -- fileUpload Button
		// rgtCol

		var ids = ["uploadBtn", "btn0", "fileUploadMessage" ];
		for ( var i = 0; i < ids.length; i++ )
		{
			if ( dijit.byId(ids[i]) )
			{
				dijit.byId(ids[i]).destroy();
			}
		}


		var handle = paneWidget.setHref('plugins/project/FileUpload/fileUpload.html');

		var fileUploadObject = this;

		handle.addOnLoad(function(){

			dojo.require("dojox.form.FileUploader");
			dojo.require("dijit.form.Button"); 
			dojo.require("dojo.parser");

			//var uploadUrl = "dojo-1.5.0/dojox/form/resources/UploadFile.php";
			var uploadUrl = Agua.cgiUrl + "upload.cgi";

			//dojo.byId("uploadedFiles").value = "";
			dojo.byId("fileToUpload").value = "";


			var button = dijit.byId('btn0');
			button.domNode.firstChild.setAttribute('class', 'fileUploadButtonSpan');
			button.domNode.firstChild.firstChild.setAttribute('class', 'fileUploadButtonSpan');
			button.domNode.firstChild.firstChild.firstChild.setAttribute('class', 'fileUploadButtonSpan');
			button.domNode.firstChild.parentNode.setAttribute('class', 'fileUploadButtonSpan');
            button.attr('iconClass', "fileUploadBrowseIcon");

			var f0 = new dojox.form.FileUploader({
				button:dijit.byId("btn0"), 
				degradable:true,
				uploadUrl:uploadUrl, 
				uploadOnChange: true, 
				selectMultipleFiles: false,
				//fileMask:fileMask,
				isDebug:true
			});



			// SET fileUploader
			fileUploadObject.fileUploader = f0;


			// OPEN LOCAL OPERATING SYSTEM FILE MANAGER
			//f0._openDialog();


			// OVERRIDE onChange METHOD IN ORDER TO
			//
			// PASS THE NAME OF THE UPLOADED FILE TO THE
			//
			// fileUp OBJECT SO THAT IT THE NAME CAN BE 
			//
			// USED IN THE fileUp.fileInfo METHOD
			f0.onChange = function(data){

				fileUploadObject.uploadedFile = data[0].name;


				dojo.forEach(data,
				function(d)
				{
					fileUploadObject.uploadedFile = d.name;
					dojo.byId("fileToUpload").value = d.name;
				});


				// summary
				//	Called after a system dialog selection has been made
				// stub to connect
				if(this.uploadOnChange) { 
					this.upload(data); 
				}else if(this.selectMultipleFiles){
					this.createFileInput();	
				}
			},


			// OVERRIDE upload FUNCTION IN ORDER T0:
			//
			// 	1. CHECK CONFIRM OVERWRITE OF EXISTING FILE
			//
			//  2. ADD path TO POST VARIABLES
			//
			f0.upload = function(data){


				var filepath = path + "/" + fileUploadObject.uploadedFile;

				// 	1. CHECK CONFIRM OVERWRITE OF EXISTING FILE
				var fileInfo = fileUploadObject.getFileInfo(filepath);

				// IF FILE EXISTS, ASK USER TO CONFIRM OVERWRITE
				var overWrite = true;
				if ( fileInfo.exists == true )
				{

					fileInfo.filepath = filepath;

					var callback = f0.executeUpload;

					overWrite = fileUploadObject.confirmOverwrite(fileInfo, callback, data);
				}

				// IF FILE DOES NOT EXIST, EXECUTE UPLOAD
				else
				{
					f0.executeUpload();

				}
			},




			f0.executeUpload = function(data)
			{

				// SHOW LOADING...
				dojo.byId("fileUploadMessage").innerHTML = "Uploading file ...";
				dojo.byId("fileUploadStatusIcon").setAttribute('class', 'fileUploadStatusIcon');


				//if ( overWrite == null || overWrite == false || ! overWrite )
				//{
				//	return;
				//}

				//  2. ADD path TO POST VARIABLES
				// DO THIS TO MAKE SURE WE DON'T KEEP ON ADDING DUPLICATE
				// HIDDEN NODES TO THE FORM
				var childNodes = f0._formNode.childNodes;
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
					f0._formNode.appendChild(userNode);

					var sessionNode = document.createElement('input');
					sessionNode.type = "hidden";
					sessionNode.name = "sessionId";
					sessionNode.value = Agua.cookie('sessionId');
					sessionNode.id = "hiddenSessionId";
					f0._formNode.appendChild(sessionNode);

					var pathNode = document.createElement('input');
					pathNode.type = "hidden";
					pathNode.name = "path";
					pathNode.value = path;
					pathNode.id = "hiddenPath";
					f0._formNode.appendChild(pathNode);
				}


				// summary
				//	Tell form to upload
				dojo.io.iframe.send({
					//url: f0.uploadUrl + "?username=admin&sessionId=1228791394.7868.158",
					url: f0.uploadUrl,
					form: f0._formNode,
					timeout: 20000,
					handleAs: "json",
					handle: dojo.hitch(f0,function(data,ioArgs,widgetRef){
						f0.onComplete(f0.selectMultipleFiles?data:[data]);								 
					})	
				});

			},


			f0.onComplete = function(dataArray){

//var t=setTimeout("javascript statement",milliseconds);

				// SHOW END LOADING
				dojo.byId("fileUploadMessage").innerHTML = "";
				dojo.byId("fileUploadStatusIcon").setAttribute('class', '');

				var filepath = fileUploadObject.path + "/" + fileUploadObject.uploadedFile;
				fileUploadObject.parentWidget.onUploadComplete(fileUploadObject.node, filepath);


				// DESTROY THIS FILE UPLOADER

				// DESTROY THE f0 OBJECT
				f0.destroyAll();

				// HIDE THE FILE UPLOAD DIALOGUE
				fileUploadObject.fileUploadDialog.hide();
			},


			dojo.connect(f0, "onProgress", function(data){
				dojo.byId("fileToUpload").value = "";
				dojo.forEach(data, function(d){
					dojo.byId("fileToUpload").value += "("+d.percent+"%) "+d.name+" \n";

					dojo.byId("fileToUpload").value = d.name;

				});
			});

			Destroy = function(){

				f0.destroyAll();
			}

		});


		handle.addOnUnload(function(){
			cleanUp();
		});


		// CREATE DIALOGUE
		var fileUploadDialog = new dojox.widget.Dialog ({
			id: dialogId,
			dimensions: [800,"auto"],
			draggable: true,
			align: 'left',
			sizeDuration: 200,
			sizeMethod:       "combine",
			viewportPadding : "125",
			refocus : false,
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
		//fileUploadDialog.show();

		// SET DIALOGUE POSITION
		//fileUploadDialog.domNode.style.left='100px';
		fileUploadDialog.domNode.style.top='100px';

		// FIXED THIS ERROR 'Error undefined running custom onLoad code'
		fileUploadDialog.dojoxDialogWrapper.appendChild(containerNode);

		containerNode.appendChild(paneNode);

		// SET fileUploadDialog
		this.fileUploadDialog = fileUploadDialog;
	},


    loadCSS : function ()
    {        

        // THIS CSS FILE PROVIDES THE ICONS AND FORMATTING FOR INDIVIDUAL FILES/DIRECTORIES
		var cssFile3 = "plugins/upload/FileUpload/fileUpload.css";
		var cssNode = document.createElement('link');
		cssNode.type = 'text/css';
		cssNode.rel = 'stylesheet';
		cssNode.href = cssFile3;
		cssNode.media = 'screen';
		cssNode.title = 'loginCSS';
		cssNode.id = "widgetStyle";
		document.getElementsByTagName("head")[0].appendChild(cssNode);
    }


}); 
