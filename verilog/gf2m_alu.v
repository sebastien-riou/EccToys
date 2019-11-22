`timescale 1ns / 1ps

`define WORD_WIDTH 256
//do not support square for now




module gf2m_alu(
    input wire clk,
    input wire reset,
    input wire add,
    input wire shl,
    input wire sto,
    input wire [`WORD_WIDTH-1:0] sbus,
    output reg [`WORD_WIDTH-1:0] dbus,
    output reg early_mz,
    output reg mz,
    output reg eq,
    output reg run,
    output reg done
    );

always @* done = sto|run;

reg add0;
reg shl0;
reg [`WORD_WIDTH-1:0] ra;
always @* begin
    add0 = add & !run;
    shl0 = shl & !run;
    ra = dbus;
end

reg [`WORD_WIDTH-1:0] next_dbus;
always @* begin
    next_dbus = ra;
    case(1'b1)
    sto : next_dbus = sbus;
    add0: next_dbus = sbus ^ ra;
    shl0: next_dbus = {sbus[`WORD_WIDTH-2:0],1'b0};
    endcase
end

always @* early_mz = sbus[`WORD_WIDTH-1];
always @* eq = ~|ra;

always @(posedge clk) begin
    if(reset) mz <= 0;
    else begin
        if(shl0|sto|add0) mz <= early_mz;
    end
end

always @(posedge clk) begin
    if(reset) dbus <= {`WORD_WIDTH{1'b0}};
    else dbus <= next_dbus;
end

always @(posedge clk) begin
    if(reset) run <= 0;
    else begin
        if(add0|shl0) run <= 1;
        if(done) run <= 0;
    end
end

endmodule
