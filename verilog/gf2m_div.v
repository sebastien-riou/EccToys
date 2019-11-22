`timescale 1ns / 1ps
`define WORD_WIDTH 256

//implemented:
//op   __-----__
//run  ___----__
//done ______-__
//out  xxxxxxxvv

module gf2m_div #(
    parameter DIV_STEP = 8
    )(
    input wire clk,
    input wire reset,
    input wire mod_div,//could also support inversion almost for free
    input wire stoy,//dividend
    input wire stox,//divider
    //input wire [`WORD_WIDTH-1:0] irreducible_poly_msb,
    input wire [`WORD_WIDTH-1:0] irreducible_poly,
    input wire [`WORD_WIDTH-1:0] sbus,
    output reg [`WORD_WIDTH-1:0] dbus,
    output reg run,
    output reg done
    );

reg run_done;
always @* done = stoy|stox|run_done;

reg [10:0] ca;
reg [10:0] cb;

reg [`WORD_WIDTH-1:0] u;
reg [`WORD_WIDTH-1:0] v;
reg [`WORD_WIDTH-1:0] a;
reg [`WORD_WIDTH-1:0] b;

always @* dbus = u;

reg [`WORD_WIDTH-1:0] next_uv;
reg [`WORD_WIDTH-1:0] next_ab;

wire [`WORD_WIDTH-1:0] one = 1;

reg a_even;
reg b_even;
reg ca_gt_cb;
reg add_both;
reg update_au;
reg update_bv;
reg end_condition;

function [`WORD_WIDTH-1:0] compute_next;
    input [`WORD_WIDTH-1:0] au;
    input [`WORD_WIDTH-1:0] bv;
    reg [`WORD_WIDTH-1:0] o;
begin
    o = 0;
    o = o ^ ((update_au | add_both) ? au : 0);
    o = o ^ ((update_bv | add_both) ? bv : 0);
    if(~o[0]) begin
        o = o>>1;
    end else begin
        o = (o^irreducible_poly)>>1;
    end
    compute_next = o;
end
endfunction

always @* begin
    a_even = ~a[0];
    b_even = ~b[0];
    ca_gt_cb = ca > cb;
    add_both = (~a_even) & (~b_even);
    update_au = (add_both &   ca_gt_cb ) | a_even;
    update_bv = (add_both & (~ca_gt_cb)) | b_even;
    end_condition = update_au & (next_ab==one);
    next_ab = compute_next(a,b);
    next_uv = compute_next(u,v);
end

always @(posedge clk,posedge reset) begin
    if(reset) begin
        ca <= 0;
        cb <= 0;
        a <= 0;
        b <= 0;
        u <= 0;
        v <= 0;
    end else begin
        if(mod_div) begin
            if(update_au) begin
                a <= next_ab;
                ca <= ca - 1;
                u <= next_uv;
            end
            if(update_bv) begin
                b <= next_ab;
                cb <= cb - 1;
                v <= next_uv;
            end
        end else begin
            if(stoy) begin
                b <= irreducible_poly;
                cb = `WORD_WIDTH;
                u <= sbus;
            end
            if(stox) begin
                a <= sbus;
                ca = `WORD_WIDTH - 1;
                v <= 0;
            end
        end
    end
end

always @(posedge clk,posedge reset) begin
    if(reset) begin
        run <= 0;
        run_done <= 0;
    end else begin
        if(mod_div) begin
            run <= 1'b1;
        end
    if(end_condition) begin
            run_done <= 1;
            run <= 1'b0;
        end
        if(run_done) begin
            run_done <= 1'b0;
        end
    end
end

endmodule
