#ifndef EXPATPP_H
#define EXPATPP_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#ifdef HAVE_EXPAT_H
#include <expat.h>
#endif

//--------------------------------------------------
// A full-blown, inlined version of Expat 2.0.1 C++ interface
//-------------------------------------------------- 
namespace expatpp {
    typedef XML_Char Char;
    typedef XML_Content Content;
    typedef XML_Content_Quant ContentQuant;
    typedef XML_Content_Type ContentType;
    typedef XML_Encoding Encoding;
    typedef XML_Error Error;
    typedef XML_Expat_Version ExpatVersion;
    typedef XML_Feature Feature;
    typedef XML_Index Index;
    typedef XML_LChar LChar;
    typedef XML_Memory_Handling_Suite MemoryHandlingSuite;
    typedef XML_ParamEntityParsing ParamEntityParsing;
    typedef XML_Parsing Parsing;
    typedef XML_ParsingStatus ParsingStatus;
    typedef XML_Size Size;
    typedef XML_Status Status;

    typedef XML_AttlistDeclHandler AttlistDeclHandler;
    typedef XML_CharacterDataHandler CharacterDataHandler;
    typedef XML_CommentHandler CommentHandler;
    typedef XML_DefaultHandler DefaultHandler;
    typedef XML_ElementDeclHandler ElementDeclHandler;
    typedef XML_EndCdataSectionHandler EndCdataSectionHandler;
    typedef XML_EndDoctypeDeclHandler EndDoctypeDeclHandler;
    typedef XML_EndElementHandler EndElementHandler;
    typedef XML_EndNamespaceDeclHandler EndNamespaceDeclHandler;
    typedef XML_EntityDeclHandler EntityDeclHandler;
    typedef XML_ExternalEntityRefHandler ExternalEntityRefHandler;
    typedef XML_NotStandaloneHandler NotStandaloneHandler;
    typedef XML_NotationDeclHandler NotationDeclHandler;
    typedef XML_ProcessingInstructionHandler ProcessingInstructionHandler;
    typedef XML_SkippedEntityHandler SkippedEntityHandler;
    typedef XML_StartCdataSectionHandler StartCdataSectionHandler;
    typedef XML_StartDoctypeDeclHandler StartDoctypeDeclHandler;
    typedef XML_StartElementHandler StartElementHandler;
    typedef XML_StartNamespaceDeclHandler StartNamespaceDeclHandler;
    typedef XML_UnknownEncodingHandler UnknownEncodingHandler;
    typedef XML_UnparsedEntityDeclHandler UnparsedEntityDeclHandler;
    typedef XML_XmlDeclHandler XmlDeclHandler;

    class Parser {
	XML_Parser parser;

    public:
	//--------------------------------------------------
	// Parser creation
	//-------------------------------------------------- 
	Parser(const Char* encoding = 0) { parser = XML_ParserCreate(encoding); } 
	Parser(const Char* encoding, Char sep) { parser = XML_ParserCreateNS(encoding, sep); }
	Parser(const Char* encoding, const MemoryHandlingSuite* ms, const Char* sep) { 
	    parser = XML_ParserCreate_MM(encoding, ms, sep); }
	Parser(Parser& e, const Char* context, const Char* encoding) {
	    parser = XML_ExternalEntityParserCreate(e.parser, context, encoding); }
	virtual ~Parser() { XML_ParserFree(parser); }
	bool reset(const Char* encoding = 0) { return XML_ParserReset(parser, encoding); }

	//--------------------------------------------------
	// Parsing
	//-------------------------------------------------- 
	int parse(const char* s, int len, int isFinal) { return XML_Parse(parser, s, len, isFinal); }
	int parseBuffer(int len, int isFinal) { return XML_ParseBuffer(parser, len, isFinal); }
	void* getBuffer(int len) { return XML_GetBuffer(parser, len); }
	int stopParser(bool resumable) { return XML_StopParser(parser, resumable); }
	int resumeParser() { return XML_ResumeParser(parser); }
	void getParsingStatus(ParsingStatus* status) { XML_GetParsingStatus(parser, status); }

	//--------------------------------------------------
	// Handler setting
	//-------------------------------------------------- 
	void setStartElementHandler(StartElementHandler start) {
	    XML_SetStartElementHandler(parser, start); }
	void setEndElementHandler(EndElementHandler end) {
	    XML_SetEndElementHandler(parser, end); }
	void setElementHandler(StartElementHandler start, EndElementHandler end) {
	    XML_SetElementHandler(parser, start, end); }
	void setCharacterDataHandler(CharacterDataHandler charhndl) {
	    XML_SetCharacterDataHandler(parser, charhndl); }
	void setProcessingInstructionHandler(ProcessingInstructionHandler proc) {
	    XML_SetProcessingInstructionHandler(parser, proc); }
	void setCommentHandler(CommentHandler cmnt) {
	    XML_SetCommentHandler(parser, cmnt); }
	void setStartCdataSectionHandler(StartCdataSectionHandler start) {
	    XML_SetStartCdataSectionHandler(parser, start); }
	void setEndCdataSectionHandler(EndCdataSectionHandler end) {
	    XML_SetEndCdataSectionHandler(parser, end); }
	void setCdataSectionHandler(StartCdataSectionHandler start, EndCdataSectionHandler end) {
	    XML_SetCdataSectionHandler(parser, start, end); }
	void setDefaultHandler(DefaultHandler hndl) {
	    XML_SetDefaultHandler(parser, hndl); }
	void setDefaultHandlerExpand(DefaultHandler hndl) {
	    XML_SetDefaultHandlerExpand(parser, hndl); }
	void setExternalEntityRefHandler(ExternalEntityRefHandler hndl) {
	    XML_SetExternalEntityRefHandler(parser, hndl); }
	void setExternalEntityRefHandlerArg(void* arg) {
	    XML_SetExternalEntityRefHandlerArg(parser, arg); }
	void setSkippedEntityHandler(SkippedEntityHandler handler) {
	    XML_SetSkippedEntityHandler(parser, handler); }
	void setUnknownEncodingHandler(UnknownEncodingHandler enchandler, void* encodingHandlerData) {
	    XML_SetUnknownEncodingHandler(parser, enchandler, encodingHandlerData); }
	void setStartNamespaceDeclHandler(StartNamespaceDeclHandler start) {
	    XML_SetStartNamespaceDeclHandler(parser, start); }
	void setEndNamespaceDeclHandler(EndNamespaceDeclHandler end) {
	    XML_SetEndNamespaceDeclHandler(parser, end); }
	void setNamespaceDeclHandler(StartNamespaceDeclHandler start, EndNamespaceDeclHandler end) {
	    XML_SetNamespaceDeclHandler(parser, start, end); }
	void setXmlDeclHandler(XmlDeclHandler xmldecl) {
	    XML_SetXmlDeclHandler(parser, xmldecl); }
	void setStartDoctypeDeclHandler(StartDoctypeDeclHandler start) {
	    XML_SetStartDoctypeDeclHandler(parser, start); }
	void setEndDoctypeDeclHandler(EndDoctypeDeclHandler end) {
	    XML_SetEndDoctypeDeclHandler(parser, end); }
	void setDoctypeDeclHandler(StartDoctypeDeclHandler start, EndDoctypeDeclHandler end) {
	    XML_SetDoctypeDeclHandler(parser, start, end); }
	void setElementDeclHandler(ElementDeclHandler eldecl) {
	    XML_SetElementDeclHandler(parser, eldecl); }
	void setAttlistDeclHandler(AttlistDeclHandler attdecl) {
	    XML_SetAttlistDeclHandler(parser, attdecl); }
	void setEntityDeclHandler(EntityDeclHandler handler) {
	    XML_SetEntityDeclHandler(parser, handler); }
	void setUnparsedEntityDeclHandler(UnparsedEntityDeclHandler h) {
	    XML_SetUnparsedEntityDeclHandler(parser, h); }
	void setNotationDeclHandler(NotationDeclHandler h) {
	    XML_SetNotationDeclHandler(parser, h); }
	void setNotStandaloneHandler(NotStandaloneHandler h) {
	    XML_SetNotStandaloneHandler(parser, h); }

	//--------------------------------------------------
	// Parse position and error reporting
	//-------------------------------------------------- 
	Error getErrorCode() { return XML_GetErrorCode(parser); }
	Index getCurrentByteIndex() { return XML_GetCurrentByteIndex(parser); }
	Size getCurrentLineNumber() { return XML_GetCurrentLineNumber(parser); }
	Size getCurrentColumnNumber() { return XML_GetCurrentColumnNumber(parser); }
	int getCurrentByteCount() { return XML_GetCurrentByteCount(parser); }
	const char* getInputContext(int* offset, int* size) {
	    return XML_GetInputContext(parser, offset, size); }

	//--------------------------------------------------
	// Miscellaneous
	//-------------------------------------------------- 
	void setUserData(void* userData) { XML_SetUserData(parser, userData); }
	void* getUserData() { return XML_GetUserData(parser); }
	void useParserAsHandlerArg() { XML_UseParserAsHandlerArg(parser); }
	int setBase(const Char* base) { return XML_SetBase(parser, base); }
	const Char* getBase() { return XML_GetBase(parser); }
	int getSpecifiedAttributeCount() { return XML_GetSpecifiedAttributeCount(parser); }
	int getIdAttributeIndex() { return XML_GetIdAttributeIndex(parser); }
	int setEncoding(const Char* encoding) { return XML_SetEncoding(parser, encoding); }
	int setParamEntityParsing(ParamEntityParsing code) { return XML_SetParamEntityParsing(parser, code); }
	Error useForeignDTD(bool useDTD) { return XML_UseForeignDTD(parser, useDTD); }
	void setReturnNSTriplet(int doNST) { XML_SetReturnNSTriplet(parser, doNST); }
	void defaultCurrent() { XML_DefaultCurrent(parser); }
	void freeContentModel(Content* model) { XML_FreeContentModel(parser, model); }
	void* memMalloc(size_t size) { return XML_MemMalloc(parser, size); }
	void* memRealloc(void* ptr, size_t size) { return XML_MemRealloc(parser, ptr, size); }
	void memFree(void* ptr) { XML_MemFree(parser, ptr); }
    };

    inline const LChar* errorString(Error code) { return XML_ErrorString(code); }
    inline const LChar* expatVersion() { return XML_ExpatVersion(); }
    inline ExpatVersion expatVersionInfo() { return XML_ExpatVersionInfo(); }
    inline const Feature* getFeatureList() { return XML_GetFeatureList(); }
}
#endif
