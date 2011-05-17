
dojo.provide("plugins.upload.FileUpload");

// DISPLAY A FILE UPLOAD DIALOGUE AND ALLOW THE USER
// TO BROWSE FILES, SELECT AND THEN UPLOAD THEM
//
// THIS MODULE (plugins.upload.FileUpload) HANDLES
// THE LOGIC AND INSTANTIATES plugins.upload.FileInput
// TO GENERATE THE 'Browse' BUTTON. THE INPUT AND
// BROWSE BUTTON NODES FROM FileInputBlind ARE SWAPPED
// FOR USER-PROVIDED NODES

// INHERITS
dojo.require("plugins.core.Common");

// DIALOGUE
dojo.require("dojox.widget.Dialog");

// FILE UPLOAD
dojo.require("plugins.upload.FileInput");	

dojo.declare( "plugins.upload.FileUpload",
	[ plugins.core.Common ],
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

	// CSS FILES
	cssFiles: [ "plugins/upload/FileUpload/fileUpload.css",
			"dojo-1.5.0/dojox/form/resources/FileInput.css" ],

	// callback FUNCTION AND DATA FROM OBJECT THAT GENERATED THE FileUpload
	onCompleteCallback : null,

	// PARENT WIDGET THAT HAS_A FileUp OBJECT
	parentWidget : null,

	// CANCEL BUTTON
	cancelButton : null,

	// REPLACE THIS NODE WITH THE 'Browse' BUTTON
	replaceNode : null,

	// USE THIS INPUT BOX AS THE FILE INPUT NODE
	inputNode : null,

	preamble: function(){
	},

 	// CONSTRUCTOR	
	constructor : function(args) {

		this.callback = args.callback;		// CALLBACK TO USE AT onComplete
		this.browseNode = args.browseNode;
		this.path = args.path;
		this.url = args.url;
		this.inputNode = args.inputNode; 		// INPUT FOR NAME OF UPLOADED FILE
		this.buttonNode = args.buttonNode; 		// CLICK THIS BUTTON TO OPEN FILE BROWSER
		this.hiddenValues = args.hiddenValues; 		// HIDDEN VALUES TO BE PASSED IN FORM


        // LOAD SORIA AND FILEPICKER CSS
        this.loadCSS();

		// CREATE FILE INPUT
		this.generateFileInput();
	},


	// GENERATE FILE UPLOAD ELEMENTS
	generateFileInput : function ()
	{

		// CREATE CONTAINER NODE TO APPEND TO TABLE ELEMENT
		//var containerNode = document.createElement('div');
		//this.buttonNode.innerHTML = '';
		//this.buttonNode.appendChild(containerNode);

		// CREATE INPUT NODE FOR fileInput
		var fileInputNode = document.createElement('input');

//return;

		// APPEND TO buttonNode IF CAN APPEND
		if ( this.buttonNode.appendChild )
		{
			this.buttonNode.appendChild(fileInputNode);
		}

		// CREATE fileInput
		var fileInput = new plugins.upload.FileInput(
			{
				url : this.url,
				label : 'Browse',
				cancelText : 'Cancel',
				blurDelay: 100,
				triggerEvent: "onchange",
				uploadMessage: "Uploading file..."
			},
			fileInputNode
		);

		// SET CALLBACK AND PATH 
		fileInput.callback = this.callback;
		fileInput.path = this.path;

		// REPLACE TITLE NODE
		fileInput.titleNode.parentNode.removeChild(fileInput.titleNode);
		fileInput.titleNode = this.buttonNode;

		// REMOVE FAKE NODE HOLDER
		fileInput.fakeNodeHolder.parentNode.removeChild(fileInput.fakeNodeHolder);

		// SET TRIGGER EVENT
		fileInput.triggerEvent = "onchange";

		//// SET BLUR LISTENER
		fileInput._blurListener = fileInput.connect(fileInput.fileInput, fileInput.triggerEvent, "_onBlur");

		// SET this.fileInput
		this.fileInput = fileInput;

		// SET CANCEL BUTTON BASE CLASS
		fileInput.set('baseClass', 'fileUpload');

		// STOP PROPAGATION WHEN CLICKED ON FILE INPUT
		dojo.connect(this.fileInput, "onclick", function (event)
			{

 				// STOP EVENT BUBBLING
				event.stopPropagation();
			}
		);

		// STOP PROPAGATION WHEN CLICK OUTSIDE FILE INPUT BOX
		dojo.connect(fileInput.domNode, "onclick", function (event)
			{

				fileInput._openedBrowser = true;

				// STOP EVENT BUBBLING
				event.stopPropagation();   
			}
		);

		// SET HIDDEN NODE VALUES
		fileInput.hiddenValues = this.hiddenValues;

		// HIDE THE OVERLAY
		fileInput.overlay.setAttribute('style', 'visibility: hidden; height: 0px');
	}

}); 
