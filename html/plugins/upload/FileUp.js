
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
	null,
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

	// CANCEL BUTTON
	cancelButton : null,

	preamble: function(){
	},

 	// CONSTRUCTOR	
	constructor : function(args) {

		this.callback = args.callback;		// CALLBACK TO USE AT onComplete
		this.path = args.path;
		this.url = args.url;
		this.node = args.buttonNode; 		// OPTION TABLE DATA NODE THE FILE IS BEING UPLOADED FOR
		this.hiddenValues = args.hiddenValues; 		// HIDDEN VALUE PAIRS TO BE ADDED TO hiddenValuesS IN FORM


		var parentNode = this.node.parentNode;

        // LOAD SORIA AND FILEPICKER CSS
        this.loadCSS();

		// CREATE FILE INPUT
		this.generateFileInput();
	},

	generateFileInput : function ()
	{
		// CREATE CONTAINER NODE TO APPEND TO TABLE ELEMENT
		var container = document.createElement('div');
		this.node.innerHTML = '';
		this.node.appendChild(container);

		// CREATE INPUT NODE FOR fileInput
		var fileInputNode = document.createElement('input');
		container.appendChild(fileInputNode);

		// CREATE fileInput
		var fileInput = new dojox.form.FileInputBlind(
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

		fileInput.triggerEvent = "onchange";
		fileInput._blurListener = fileInput.connect(fileInput.fileInput, fileInput.triggerEvent, "_onBlur");

		// SET AS this.fileInput
		this.fileInput = fileInput;


		// SET CANCEL BUTTON BASE CLASS
		fileInput.attr('baseClass', 'fileUpload');

		// AFTER MOUSEOVER
		//<span class="dijit dijitReset dijitLeft dijitInline adminApplicationTableRowcancel fileUploadCancelButton fileUploadCancelButtonHover dijitHover" dojoattachevent="ondijitclick:_onButtonClick,onmouseenter:_onMouse,onmouseleave:_onMouse,onmousedown:_onMouse" style="width: auto;" widgetid="dojoUnique1">
		//
		// BEFORE MOUSEOVER 
		//<span class="dijit dijitReset dijitLeft dijitInline dijitButton adminApplicationTableRowcancel" dojoattachevent="ondijitclick:_onButtonClick,onmouseenter:_onMouse,onmouseleave:_onMouse,onmousedown:_onMouse" style="width: auto;" widgetid="dojoUnique1">
		//
		// WITH focus()
		//<span class="dijit dijitReset dijitLeft dijitInline adminApplicationTableRowcancel fileUploadCancelButton" dojoattachevent="ondijitclick:_onButtonClick,onmouseenter:_onMouse,onmouseleave:_onMouse,onmousedown:_onMouse" style="width: auto;" widgetid="dojoUnique1">

		// FOCUS TO APPLY CSS STYLING
		fileInput.focusNode.focus();

		// STOP PROPAGATION WHEN CLICK OUTSIDE FILE INPUT BOX
		dojo.connect(fileInput.fileInput, "onclick", function (event)
			{
				// STOP EVENT BUBBLING
				event.stopPropagation();   
			}
		);

		// STOP PROPAGATION WHEN CLICK OUTSIDE FILE INPUT BOX
		dojo.connect(fileInput.domNode, "onclick", function (event)
			{
				// STOP EVENT BUBBLING
				event.stopPropagation();   
			}
		);

		var blindInput = fileInput.fakeNodeHolder.childNodes[1];
		console.dir(blindInput);




		// REMOVE THE 'CANCEL' NODE
		fileInput.fakeNodeHolder.removeChild(fileInput.cancelNode);






		dojo.connect(document, "onkeypress", function(event){

			var key = event.charOrCode;

			// SAVE IF 'RETURN' KEY PRESSED
			if ( key == 13 )
			{
				//workflowObject.saveOptionValue(this);
				//workflowObject.revertOptionValue(this);
				//return;
			}
		});


		// SET RETURN KEY LISTENER
		//var workflowObject = this;
		dojo.connect(blindInput, "onkeypress", function(event){

			var key = event.charOrCode;

			// SAVE IF 'RETURN' KEY PRESSED
			if ( key == 13 )
			{
				//workflowObject.saveOptionValue(this);
				//workflowObject.revertOptionValue(this);
				//return;
			}

			// RETURN IF 'ESCAPE' KEY IS PRESSED
			if (key == dojo.keys.ESCAPE)
			{
				// STOP EVENT BUBBLING
				event.stopPropagation();   

				// DO onComplete TO QUIT UPLOAD
				fileUpObject.fileInput.onComplete();
			}

			// DO NOTHING IF 'BACKSPACE' KEY IS PRESSED
			if (key == dojo.keys.BACKSPACE)
			{
				//return;
			}
		});



		// CANCEL BUTTON ONCLICK
		var fileUpObject = this;
		dojo.connect(fileInput.cancelNode, "onclick", function(event)
			{

				// STOP EVENT BUBBLING
				event.stopPropagation();   

				// DO onComplete TO QUIT UPLOAD
				fileUpObject.fileInput.onComplete();
			}
		);



		// SET HIDDEN NODE VALUES
		fileInput.hiddenValues = this.hiddenValues;

		// DO NOT START UP FILE INPUT (WILL MAKE TEXT INPUT DISAPPEAR)
		////////fileInput.startup();

		//HIDE THE PROGRESS OVERLAY FOR NOW
		//<div class="dijitFileInput " widgetid="dojox_form_FileInputBlind_0">
		//	<input id="dojox_form_FileInputBlind_0" class="dijitFileInputReal" type="file" dojoattachpoint="fileInput" name="uploadFile" style="left: -142px;"/>
		//		<div class="dijitFakeInput" dojoattachpoint="fakeNodeHolder">
		//			<div class="dijitInline dijitFileInputText" dojoattachpoint="titleNode"/>
		//			<div class="dijitInline dijitFileInputButton" dojoattachevent="onclick:reset" dojoattachpoint="cancelNode"/>
		//		</div>
		//	<div class="dijitProgressOverlay" dojoattachpoint="overlay"> </div>
		//</div>

		fileInput.overlay.setAttribute('style', 'visibility: hidden; height: 0px');

		// OVERRIDE OF SENDFILE TO ADD FILE PATH INFO
		// summary: triggers the chain of events needed to upload a file in the background.
		fileInput._sendFile = function(e){

			// IF THE INPUT ELEMENT HAS NO VALUE, ABORT THE SEND AND DESTROY THIS OBJECT

			if (this._sent || this._sending || !this.fileInput.value)
			{
				return;
			}

			this.overlay.setAttribute('style', 'visibility: visible; height: 20px; width: auto;');

			this._sending = true;

			dojo.style(this.fakeNodeHolder,"display","none");
			dojo.style(this.overlay,{
				opacity:0,
				display:"block"
			});

			this.setMessage(this.uploadMessage);
			dojo.fadeIn({ node: this.overlay, duration:this.duration }).play();

			var _newForm; 
			if(dojo.isIE){
				// just to reiterate, IE is a steaming pile of code. 
				_newForm = document.createElement('<form enctype="multipart/form-data" method="post">');
				_newForm.encoding = "multipart/form-data";	
			}
			else
			{
				// this is how all other sane browsers do it
				_newForm = document.createElement('form');
				_newForm.setAttribute("enctype","multipart/form-data");
			}

			_newForm.appendChild(this.fileInput);
			dojo.body().appendChild(_newForm);


			// LATER: FINISH THIS - AUTOMATE CREATION OF HIDDEN NODES
			var hiddenValues = this.hiddenValues;
			for ( var i = 0; i < hiddenValues.length; i++ )
			{


			}


			var userNode = document.createElement('input');
			userNode.type = "hidden";
			userNode.name = "username";
			userNode.value = Agua.cookie('username');
			userNode.id = "hiddenUsername";
			_newForm.appendChild(userNode);

			var sessionNode = document.createElement('input');
			sessionNode.type = "hidden";
			sessionNode.name = "sessionId";
			sessionNode.value = Agua.cookie('sessionId');
			sessionNode.id = "hiddenSessionId";
			_newForm.appendChild(sessionNode);

			var pathNode = document.createElement('input');
			pathNode.type = "hidden";
			pathNode.name = "path";
			pathNode.value = fileUpObject.path;
			pathNode.id = "hiddenPath";
			_newForm.appendChild(pathNode);



			dojo.io.iframe.send({
				url: this.url,
				form: _newForm,
				timeout : 1000,
				handleAs: "json",
				handle: dojo.hitch(this,"onComplete")
			});
		};


		fileInput.onComplete = function(){


			// REMOVE 'Uploading file' MESSAGE
			fileInput.overlay.setAttribute('style', 'visibility: hidden;');

			// SET THE UPLOADED FILE NAME IN THE TABLE ELEMENT 

			if ( fileInput.fileInput.value == null || fileInput.fileInput.value == '' )
			{
				fileUpObject.callback();
			}
			else
			{
				var filepath = fileUpObject.path + "/" + fileInput.fileInput.value;

				// DO CALLBACK TO UPDATE THE WORKFLOW AND THE OPTIONS TABLE
				fileUpObject.callback(filepath);				
			}

			// DESTROY THE fileInput OBJECT

			this.destroy();
		};
	},


    loadCSS : function ()
    {        

        // ADDITIONAL STYLE FORMATTING FOR fileUp CLASS
		var cssFile3 = "plugins/upload/FileUpload/fileUpload.css";
		var cssNode = document.createElement('link');
		cssNode.type = 'text/css';
		cssNode.rel = 'stylesheet';
		cssNode.href = cssFile3;
		cssNode.media = 'screen';
		cssNode.title = 'fileUpload';
		cssNode.id = "widgetStyle";
		document.getElementsByTagName("head")[0].appendChild(cssNode);


		// ORIGINAL CSS FOR INHERITED FileInput CLASS
		var cssFile3 = "dojo-1.5.0/dojox/form/resources/FileInput.css";
		var cssNode = document.createElement('link');
		cssNode.type = 'text/css';
		cssNode.rel = 'stylesheet';
		cssNode.href = cssFile3;
		cssNode.media = 'screen';
		cssNode.title = 'fileInput';
		cssNode.id = "widgetStyle";
		document.getElementsByTagName("head")[0].appendChild(cssNode);
    }

}); 

